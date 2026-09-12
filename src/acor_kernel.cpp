#include "acor_kernel.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Independent diagnostic prototype. Reads arrive one at a time as strict TSV:
// dataset, cluster, read_id, sequence, is_first, is_last, total_reads, design_length.
// Candidate A light implementation. The stopping rule and its exact cumulative
// residual evidence are unchanged. Engineering changes are limited to an
// allocation-free ED kernel, memoized per-hypothesis updates, and compact
// cluster-level diagnostics instead of a sequence-heavy per-prefix trace.

enum class ConstraintMode {
    None,
    HardEndpoint,
    SoftLength,
};

enum class OptimizationMode {
    Full,
    Banded,
    FinalOnlyExact,
    SparseIncrementalExactFinal,
};

struct BandPlan {
    int start_transition = 0;
    int end_transition = 0;
    int center_transition = 0;
};

struct Config {
    std::string name;
    int k = 7;
    int beam_width = 8;
    std::string length_prior = "median_band";
};

static constexpr uint32_t NO_PATH = std::numeric_limits<uint32_t>::max();

struct PathNode {
    uint32_t parent = NO_PATH;
    uint32_t start_code = 0;
    uint32_t context = 0;
    uint32_t length = 0;
    int transitions = 0;
    uint8_t appended_base = 0;
    std::array<uint32_t, 4> children{NO_PATH, NO_PATH, NO_PATH, NO_PATH};
};

struct BeamState {
    uint32_t path = NO_PATH;
    double cumulative_score = -std::numeric_limits<double>::infinity();
};

struct Candidate {
    std::string sequence;
    double score = -std::numeric_limits<double>::infinity();
};

struct SparseCandidateState {
    std::string sequence;
    uint32_t start_code = 0;
    uint32_t end_code = 0;
    int transitions = 0;
    double transition_sum = 0.0;
    double score = -std::numeric_limits<double>::infinity();
    std::unordered_map<uint32_t, std::array<uint32_t, 4>> edge_occurrences;
    std::unordered_map<uint32_t, std::vector<uint32_t>> context_positions;
    std::unordered_map<uint32_t, double> edge_score_cache;
};

struct PathCandidate {
    uint32_t path = NO_PATH;
    double score = -std::numeric_limits<double>::infinity();
};

#ifdef ACOR_FRONTIER_AUDIT
static std::ostream* frontier_audit_output = nullptr;
static uint64_t frontier_audit_event = 0;

static void set_frontier_audit_output(std::ostream* output) {
    frontier_audit_output = output;
    frontier_audit_event = 0;
    if (frontier_audit_output != nullptr) {
        *frontier_audit_output
            << "event\treads_arrived\tlabel\trank\tsequence\tscore_decimal\tscore_hex"
               "\ttie_sequence\torder_summary\n";
    }
}

static uint64_t frontier_score_bits(double score) {
    uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(score));
    std::memcpy(&bits, &score, sizeof(bits));
    return bits;
}
#endif

struct ConsensusResult {
    std::vector<Candidate> ranked;
    double median_read_length = 0.0;
    int minimum_length = 0;
    int maximum_length = 0;
    std::string selection_mode = "unset";
    int selected_length_delta = 0;
    bool fallback_to_arm0 = false;
    int exact_m_candidate_count = 0;
    int soft_window_candidate_count = 0;
    double final_decode_seconds = 0.0;
    int band_start_transition = 0;
    int band_end_transition = 0;
    int band_center_transition = 0;
    int dirty_context_count = 0;
    int sparse_candidates_refreshed = 0;
    int sparse_dependencies_updated = 0;
    int sparse_repairs_generated = 0;
};

static char canonical_base(char ch) {
    switch (ch) {
        case 'A': case 'a': return 'A';
        case 'C': case 'c': return 'C';
        case 'G': case 'g': return 'G';
        case 'T': case 't': return 'T';
        default: return 'N';
    }
}

static int base_code(char ch) {
    switch (ch) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

static char code_base(int code) {
    static constexpr char BASES[] = {'A', 'C', 'G', 'T'};
    return BASES[code & 3];
}

static std::string canonical(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (char ch : sequence) result.push_back(canonical_base(ch));
    return result;
}

static int edit_distance_reference(const std::string& lhs, const std::string& rhs) {
    std::vector<int> previous(rhs.size() + 1U);
    std::vector<int> current(rhs.size() + 1U);
    std::iota(previous.begin(), previous.end(), 0);
    for (size_t i = 1; i <= lhs.size(); ++i) {
        current[0] = static_cast<int>(i);
        for (size_t j = 1; j <= rhs.size(); ++j) {
            current[j] = std::min({
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + (lhs[i - 1U] != rhs[j - 1U]),
            });
        }
        previous.swap(current);
    }
    return previous.back();
}

static int myers_64(const char* pattern, size_t pattern_size, const char* text, size_t text_size) {
    std::array<uint64_t, 256> equality{};
    for (size_t i = 0; i < pattern_size; ++i) {
        equality[static_cast<unsigned char>(pattern[i])] |= uint64_t{1} << i;
    }
    uint64_t positive = ~uint64_t{0};
    uint64_t negative = 0;
    int score = static_cast<int>(pattern_size);
    const uint64_t highest = uint64_t{1} << (pattern_size - 1U);
    for (size_t i = 0; i < text_size; ++i) {
        const uint64_t equal = equality[static_cast<unsigned char>(text[i])];
        const uint64_t vertical = equal | negative;
        const uint64_t horizontal = (((equal & positive) + positive) ^ positive) | equal;
        uint64_t positive_horizontal = negative | ~(horizontal | positive);
        uint64_t negative_horizontal = positive & horizontal;
        if (positive_horizontal & highest) ++score;
        if (negative_horizontal & highest) --score;
        positive_horizontal = (positive_horizontal << 1U) | 1U;
        negative_horizontal <<= 1U;
        positive = negative_horizontal | ~(vertical | positive_horizontal);
        negative = positive_horizontal & vertical;
    }
    return score;
}

static int myers_blocks(const char* pattern, size_t pattern_size, const char* text, size_t text_size) {
    static constexpr size_t MAX_BLOCKS = 16U;
    const size_t blocks = (pattern_size + 63U) / 64U;
    if (blocks == 0U || blocks > MAX_BLOCKS) return -1;
    std::array<std::array<uint64_t, MAX_BLOCKS>, 5> equality;
    auto symbol = [](char base) -> size_t {
        switch (base) {
            case 'A': return 0U;
            case 'C': return 1U;
            case 'G': return 2U;
            case 'T': return 3U;
            default: return 4U;
        }
    };
    for (auto& symbol_blocks : equality) {
        for (size_t block = 0; block < blocks; ++block) symbol_blocks[block] = 0U;
    }
    for (size_t i = 0; i < pattern_size; ++i) {
        equality[symbol(pattern[i])][i / 64U] |= uint64_t{1} << (i % 64U);
    }
    std::array<uint64_t, MAX_BLOCKS> positive;
    std::array<uint64_t, MAX_BLOCKS> negative;
    for (size_t block = 0; block < blocks; ++block) {
        positive[block] = ~uint64_t{0};
        negative[block] = 0U;
    }
    int score = static_cast<int>(pattern_size);
    const uint64_t highest = uint64_t{1} << ((pattern_size - 1U) % 64U);
    for (size_t text_index = 0; text_index < text_size; ++text_index) {
        uint64_t addition_carry = 0U;
        uint64_t positive_shift_carry = 1U;
        uint64_t negative_shift_carry = 0U;
        for (size_t block = 0; block < blocks; ++block) {
            const uint64_t equal = equality[symbol(text[text_index])][block];
            const uint64_t vertical = equal | negative[block];
            const unsigned __int128 sum =
                static_cast<unsigned __int128>(equal & positive[block]) +
                positive[block] + addition_carry;
            addition_carry = static_cast<uint64_t>(sum >> 64U);
            const uint64_t horizontal =
                (static_cast<uint64_t>(sum) ^ positive[block]) | equal;
            uint64_t positive_horizontal = negative[block] | ~(horizontal | positive[block]);
            uint64_t negative_horizontal = positive[block] & horizontal;
            if (block + 1U == blocks) {
                if (positive_horizontal & highest) ++score;
                if (negative_horizontal & highest) --score;
            }
            const uint64_t next_positive_carry = positive_horizontal >> 63U;
            const uint64_t next_negative_carry = negative_horizontal >> 63U;
            positive_horizontal = (positive_horizontal << 1U) | positive_shift_carry;
            negative_horizontal = (negative_horizontal << 1U) | negative_shift_carry;
            positive_shift_carry = next_positive_carry;
            negative_shift_carry = next_negative_carry;
            positive[block] = negative_horizontal | ~(vertical | positive_horizontal);
            negative[block] = positive_horizontal & vertical;
        }
    }
    return score;
}

static int shifted_alignment_upper(
    const std::string& lhs, size_t lhs_begin, size_t lhs_size,
    const std::string& rhs, size_t rhs_begin, size_t rhs_size) {
    int best = static_cast<int>(std::max(lhs_size, rhs_size));
    const int length_delta = static_cast<int>(lhs_size) - static_cast<int>(rhs_size);
    std::array<int, 18> shifts{};
    size_t count = 0;
    for (int delta = -4; delta <= 4; ++delta) shifts[count++] = delta;
    for (int delta = -4; delta <= 4; ++delta) shifts[count++] = length_delta + delta;
    for (size_t shift_index = 0; shift_index < count; ++shift_index) {
        const int shift = shifts[shift_index];
        const size_t left_skip = shift > 0 ? static_cast<size_t>(shift) : 0U;
        const size_t right_skip = shift < 0 ? static_cast<size_t>(-shift) : 0U;
        if (left_skip > lhs_size || right_skip > rhs_size) continue;
        const size_t overlap = std::min(lhs_size - left_skip, rhs_size - right_skip);
        int cost = static_cast<int>(left_skip + right_skip);
        for (size_t i = 0; i < overlap; ++i) {
            cost += lhs[lhs_begin + left_skip + i] != rhs[rhs_begin + right_skip + i];
        }
        cost += static_cast<int>((lhs_size - left_skip - overlap) + (rhs_size - right_skip - overlap));
        best = std::min(best, cost);
    }
    return best;
}

static int edit_distance(const std::string& lhs, const std::string& rhs) {
    size_t begin = 0;
    while (begin < lhs.size() && begin < rhs.size() && lhs[begin] == rhs[begin]) ++begin;
    size_t lhs_end = lhs.size();
    size_t rhs_end = rhs.size();
    while (lhs_end > begin && rhs_end > begin && lhs[lhs_end - 1U] == rhs[rhs_end - 1U]) {
        --lhs_end;
        --rhs_end;
    }
    const size_t lhs_size = lhs_end - begin;
    const size_t rhs_size = rhs_end - begin;
    if (lhs_size == 0U) return static_cast<int>(rhs_size);
    if (rhs_size == 0U) return static_cast<int>(lhs_size);
    const bool lhs_is_pattern = lhs_size <= rhs_size;
    const char* pattern = lhs_is_pattern ? lhs.data() + begin : rhs.data() + begin;
    const char* text = lhs_is_pattern ? rhs.data() + begin : lhs.data() + begin;
    const size_t pattern_size = lhs_is_pattern ? lhs_size : rhs_size;
    const size_t text_size = lhs_is_pattern ? rhs_size : lhs_size;
    if (pattern_size <= 64U) {
        return myers_64(pattern, pattern_size, text, text_size);
    }
    const int bit_parallel = myers_blocks(pattern, pattern_size, text, text_size);
    if (bit_parallel >= 0) {
        return bit_parallel;
    }

    // A concrete edit script supplies an exact upper bound. Ukkonen's diagonal
    // band at that width cannot exclude any path with cost <= the bound.
    const int band = shifted_alignment_upper(lhs, begin, lhs_size, rhs, begin, rhs_size);
    const int infinity = band + static_cast<int>(lhs_size + rhs_size) + 1;
    static thread_local std::vector<int> previous;
    static thread_local std::vector<int> current;
    previous.resize(rhs_size + 1U);
    current.resize(rhs_size + 1U);
    for (size_t j = 0; j <= rhs_size; ++j) {
        previous[j] = j <= static_cast<size_t>(band) ? static_cast<int>(j) : infinity;
    }
    for (size_t i = 1; i <= lhs_size; ++i) {
        const size_t start = i > static_cast<size_t>(band) ? i - static_cast<size_t>(band) : 1U;
        const size_t end = std::min(rhs_size, i + static_cast<size_t>(band));
        current[0] = i <= static_cast<size_t>(band) ? static_cast<int>(i) : infinity;
        if (start > 1U) current[start - 1U] = infinity;
        for (size_t j = start; j <= end; ++j) {
            current[j] = std::min({
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + (lhs[begin + i - 1U] != rhs[begin + j - 1U]),
            });
        }
        if (end < rhs_size) current[end + 1U] = infinity;
        previous.swap(current);
    }
    return previous.back();
}

static void edit_distance_self_test() {
    uint64_t state = 0x6a09e667f3bcc909ULL;
    auto random_u32 = [&]() -> uint32_t {
        state ^= state << 13U;
        state ^= state >> 7U;
        state ^= state << 17U;
        return static_cast<uint32_t>(state >> 16U);
    };
    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    for (int trial = 0; trial < 256; ++trial) {
        const size_t lhs_size = static_cast<size_t>(random_u32() % 181U);
        const size_t rhs_size = static_cast<size_t>(random_u32() % 181U);
        std::string lhs(lhs_size, 'A');
        std::string rhs(rhs_size, 'A');
        for (char& base : lhs) base = bases[random_u32() & 3U];
        for (char& base : rhs) base = bases[random_u32() & 3U];
        const int expected = edit_distance_reference(lhs, rhs);
        const int observed = edit_distance(lhs, rhs);
        if (expected != observed) {
            throw std::runtime_error("optimized edit-distance self-test failed");
        }
    }
}

struct ResidualCostCache {
    size_t reads_scored = 0;
    int64_t cumulative_ed = 0;
};

static uint64_t fnv1a64(const std::string& value) {
    uint64_t hash = 1469598103934665603ULL;
    for (unsigned char byte : value) {
        hash ^= static_cast<uint64_t>(byte);
        hash *= 1099511628211ULL;
    }
    return hash;
}

static std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

static Config parse_config(const std::string& name) {
    Config config;
    config.name = name;
    const std::string prefix = "kmc_k";
    if (name.rfind(prefix, 0) != 0) throw std::runtime_error("invalid config prefix: " + name);
    size_t beam_marker = name.find("_b", prefix.size());
    if (beam_marker == std::string::npos) throw std::runtime_error("missing beam width: " + name);
    size_t prior_marker = name.find('_', beam_marker + 2);
    if (prior_marker == std::string::npos) throw std::runtime_error("missing length prior: " + name);
    auto parse_component = [&](size_t begin, size_t end, const char* label) -> int {
        const std::string_view text(name.data() + begin, end - begin);
        int value = 0;
        const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
        if (text.empty() || parsed.ec != std::errc() || parsed.ptr != text.data() + text.size()) {
            throw std::runtime_error(std::string("invalid ") + label + ": " + name);
        }
        return value;
    };
    config.k = parse_component(prefix.size(), beam_marker, "k");
    config.beam_width = parse_component(beam_marker + 2, prior_marker, "beam width");
    config.length_prior = name.substr(prior_marker + 1);
    if ((config.k != 5 && config.k != 7 && config.k != 9) ||
        (config.beam_width != 4 && config.beam_width != 8 && config.beam_width != 16) ||
        (config.length_prior != "median_band" && config.length_prior != "lognormal")) {
        throw std::runtime_error("unsupported config: " + name);
    }
    return config;
}

static double median(std::vector<int> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    if (n % 2) return static_cast<double>(values[n / 2]);
    return 0.5 * static_cast<double>(values[n / 2 - 1] + values[n / 2]);
}

static int compare_finite_scores_descending(double lhs, double rhs) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        throw std::runtime_error("non-finite score in ordering comparator");
    }
    constexpr long double SCALE = 1000000000000.0L;
    const long double lhs_scaled = static_cast<long double>(lhs) * SCALE;
    const long double rhs_scaled = static_cast<long double>(rhs) * SCALE;
    if (lhs_scaled < static_cast<long double>(std::numeric_limits<int64_t>::min()) ||
        lhs_scaled > static_cast<long double>(std::numeric_limits<int64_t>::max()) ||
        rhs_scaled < static_cast<long double>(std::numeric_limits<int64_t>::min()) ||
        rhs_scaled > static_cast<long double>(std::numeric_limits<int64_t>::max())) {
        throw std::runtime_error("score outside comparator quantization range");
    }
    const int64_t lhs_bucket = static_cast<int64_t>(std::llround(lhs_scaled));
    const int64_t rhs_bucket = static_cast<int64_t>(std::llround(rhs_scaled));
    if (lhs_bucket > rhs_bucket) return -1;
    if (lhs_bucket < rhs_bucket) return 1;
    return 0;
}

static bool candidate_better(const Candidate& lhs, const Candidate& rhs) {
    const int score_order = compare_finite_scores_descending(lhs.score, rhs.score);
    if (score_order != 0) return score_order < 0;
    return lhs.sequence < rhs.sequence;
}

static bool path_lex_less_equal_length(
    const std::vector<PathNode>& arena, uint32_t lhs, uint32_t rhs) {
    if (lhs == rhs) return false;
    uint32_t left = lhs;
    uint32_t right = rhs;
    while (arena[left].parent != arena[right].parent) {
        left = arena[left].parent;
        right = arena[right].parent;
    }
    if (arena[left].parent == NO_PATH) return arena[left].start_code < arena[right].start_code;
    return arena[left].appended_base < arena[right].appended_base;
}

static std::string materialize_path(
    const std::vector<PathNode>& arena, uint32_t path, int k) {
    std::string sequence(static_cast<size_t>(arena[path].length), 'A');
    uint32_t current = path;
    size_t index = sequence.size();
    while (arena[current].parent != NO_PATH) {
        sequence[--index] = code_base(static_cast<int>(arena[current].appended_base));
        current = arena[current].parent;
    }
    uint32_t code = arena[current].start_code;
    for (int offset = k - 1; offset >= 0; --offset) {
        sequence[static_cast<size_t>(offset)] = code_base(static_cast<int>(code & 3U));
        code >>= 2U;
    }
    return sequence;
}

class DenseCounter {
public:
    explicit DenseCounter(size_t size = 0) : values_(size, 0U) {}

    void resize(size_t size) {
        values_.assign(size, 0U);
        touched_.clear();
    }

    void clear() {
        for (uint32_t index : touched_) values_[index] = 0U;
        touched_.clear();
    }

    void increment(uint32_t index) {
        if (values_[index] == 0U) touched_.push_back(index);
        if (values_[index] == std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("dense counter overflow");
        }
        ++values_[index];
    }

    uint32_t get(uint32_t index) const { return values_[index]; }
    bool empty() const { return touched_.empty(); }
    size_t size() const { return touched_.size(); }
    const std::vector<uint32_t>& touched() const { return touched_; }

private:
    std::vector<uint32_t> values_;
    std::vector<uint32_t> touched_;
};

class IncrementalKmcBeam {
public:
    IncrementalKmcBeam(
        Config config,
        ConstraintMode constraint_mode,
        OptimizationMode optimization_mode,
        int band_width)
        : config_(std::move(config)),
          constraint_mode_(constraint_mode),
          optimization_mode_(optimization_mode),
          band_width_(band_width) {
        if (optimization_mode_ == OptimizationMode::Banded && band_width_ <= 0) {
            throw std::runtime_error("banded mode requires positive band width");
        }
        mask_ = (static_cast<uint32_t>(1) << (2 * config_.k)) - 1U;
        const size_t node_slots = static_cast<size_t>(1U) << (2 * config_.k);
        node_counts_.resize(node_slots);
        outgoing_totals_.resize(node_slots);
        start_counts_.resize(node_slots);
        end_counts_.resize(node_slots);
        edge_counts_.resize(node_slots * 4U);
        transition_score_cache_.assign(
            node_slots * 4U, -std::numeric_limits<double>::infinity());
        transition_score_cache_epoch_.assign(node_slots, 0U);
        dirty_context_epoch_.assign(node_slots, 0U);
    }

    void reset(int design_length) {
        if (design_length < config_.k || design_length > std::numeric_limits<int>::max() - 5) {
            throw std::runtime_error("design_length outside safe kernel range");
        }
        node_counts_.clear();
        edge_counts_.clear();
        outgoing_totals_.clear();
        start_counts_.clear();
        end_counts_.clear();
        read_counts_.clear();
        lower_lengths_ = {};
        upper_lengths_ = {};
        log_length_sum_ = 0.0;
        log_length_sq_sum_ = 0.0;
        previous_candidate_.clear();
        candidate_switches_ = 0;
        reads_arrived_ = 0;
        design_length_ = design_length;
        exact_m_read_count_ = 0;
        graph_update_seconds_ = 0.0;
        beam_decode_seconds_ = 0.0;
        start_paths_.clear();
        path_arena_.clear();
        path_nodes_created_ = 0;
        path_cache_hits_ = 0;
        sparse_pool_.clear();
        dirty_contexts_.clear();
        sparse_candidates_refreshed_total_ = 0;
        sparse_dependencies_updated_total_ = 0;
        sparse_repairs_generated_total_ = 0;
        total_start_observations_ = 0;
    }

    ConsensusResult add_read(const std::string& sequence, bool is_last) {
        if (sequence.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("read length exceeds integer range");
        }
        if (reads_arrived_ == std::numeric_limits<size_t>::max()) {
            throw std::runtime_error("read arrival counter overflow");
        }
        ++reads_arrived_;
        add_length(static_cast<int>(sequence.size()));
        const double log_length = std::log(static_cast<double>(std::max<size_t>(1, sequence.size())));
        log_length_sum_ += log_length;
        log_length_sq_sum_ += log_length * log_length;
        uint32_t& sequence_count = read_counts_[sequence];
        if (sequence_count == std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("read multiplicity overflow");
        }
        ++sequence_count;
        if (static_cast<int>(sequence.size()) == design_length_) ++exact_m_read_count_;
        const auto graph_started = std::chrono::steady_clock::now();
        update_graph(sequence);
        graph_update_seconds_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - graph_started).count();
        if (optimization_mode_ == OptimizationMode::FinalOnlyExact && !is_last) {
            ConsensusResult deferred;
            deferred.median_read_length = current_median();
            const auto bounds = length_bounds(deferred.median_read_length);
            deferred.minimum_length = bounds.first;
            deferred.maximum_length = bounds.second;
            deferred.selection_mode = "deferred_exact_final_decode";
            return deferred;
        }
        const auto decode_started = std::chrono::steady_clock::now();
        ConsensusResult result;
        if (optimization_mode_ == OptimizationMode::SparseIncrementalExactFinal && !is_last) {
            result = recompute_sparse(sequence);
        } else {
            result = constraint_mode_ == ConstraintMode::SoftLength && is_last
                ? recompute_soft() : recompute();
        }
        const double decode_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - decode_started).count();
        beam_decode_seconds_ += decode_seconds;
        if (is_last) {
            result.final_decode_seconds = decode_seconds;
        }
        const std::string current = result.ranked.empty() ? std::string() : result.ranked.front().sequence;
        if (!previous_candidate_.empty() && current != previous_candidate_) ++candidate_switches_;
        previous_candidate_ = current;
        return result;
    }

    size_t node_count() const { return node_counts_.size(); }
    size_t edge_count() const { return edge_counts_.size(); }
    size_t reads_arrived() const { return reads_arrived_; }
    size_t candidate_switches() const { return candidate_switches_; }
    const std::string& current_candidate() const { return previous_candidate_; }
    int design_length() const { return design_length_; }
    size_t exact_m_read_count() const { return exact_m_read_count_; }
    ConsensusResult soft_decode_now() { return recompute_soft(); }
    double graph_update_seconds() const { return graph_update_seconds_; }
    double beam_decode_seconds() const { return beam_decode_seconds_; }
    uint64_t path_nodes_created() const { return path_nodes_created_; }
    uint64_t path_cache_hits() const { return path_cache_hits_; }
    uint64_t sparse_candidates_refreshed_total() const { return sparse_candidates_refreshed_total_; }
    uint64_t sparse_dependencies_updated_total() const { return sparse_dependencies_updated_total_; }
    uint64_t sparse_repairs_generated_total() const { return sparse_repairs_generated_total_; }
    ConstraintMode constraint_mode() const { return constraint_mode_; }
    OptimizationMode optimization_mode() const { return optimization_mode_; }
    int band_width() const { return band_width_; }
    bool emit_dense_diagnostics() const {
        return optimization_mode_ != OptimizationMode::FinalOnlyExact;
    }
#ifdef ACOR_TESTING
    bool test_epoch_wrap_reset() {
        current_update_epoch_ = std::numeric_limits<uint32_t>::max();
        dirty_context_epoch_[0] = 17U;
        transition_score_cache_epoch_[0] = 17U;
        transition_score_cache_[0] = 123.0;
        begin_dirty_epoch();
        return current_update_epoch_ == 1U
            && dirty_context_epoch_[0] == 0U
            && transition_score_cache_epoch_[0] == std::numeric_limits<uint32_t>::max()
            && std::isinf(transition_score_cache_[0])
            && transition_score_cache_[0] < 0.0;
    }
#endif

private:
    Config config_;
    ConstraintMode constraint_mode_ = ConstraintMode::None;
    OptimizationMode optimization_mode_ = OptimizationMode::Full;
    int band_width_ = 0;
    int design_length_ = 0;
    uint32_t mask_ = 0;
    DenseCounter node_counts_;
    DenseCounter edge_counts_;
    DenseCounter outgoing_totals_;
    DenseCounter start_counts_;
    DenseCounter end_counts_;
    mutable std::vector<double> transition_score_cache_;
    mutable std::vector<uint32_t> transition_score_cache_epoch_;
    std::unordered_map<std::string, uint32_t> read_counts_;
    std::priority_queue<int> lower_lengths_;
    std::priority_queue<int, std::vector<int>, std::greater<int>> upper_lengths_;
    double log_length_sum_ = 0.0;
    double log_length_sq_sum_ = 0.0;
    std::string previous_candidate_;
    size_t candidate_switches_ = 0;
    size_t reads_arrived_ = 0;
    size_t exact_m_read_count_ = 0;
    double graph_update_seconds_ = 0.0;
    double beam_decode_seconds_ = 0.0;
    std::unordered_map<uint32_t, uint32_t> start_paths_;
    std::vector<PathNode> path_arena_;
    uint64_t path_nodes_created_ = 0;
    uint64_t path_cache_hits_ = 0;
    std::vector<SparseCandidateState> sparse_pool_;
    std::vector<uint32_t> dirty_contexts_;
    std::vector<uint32_t> dirty_context_epoch_;
    uint32_t current_update_epoch_ = 0;
    uint64_t sparse_candidates_refreshed_total_ = 0;
    uint64_t sparse_dependencies_updated_total_ = 0;
    uint64_t sparse_repairs_generated_total_ = 0;
    uint64_t total_start_observations_ = 0;

#ifdef ACOR_FRONTIER_AUDIT
    void audit_frontier_rows(
        const char* label,
        const std::vector<std::pair<std::string, double>>& rows) const {
        if (frontier_audit_output == nullptr) return;
        const uint64_t event = ++frontier_audit_event;
        std::string summary_material;
        for (const auto& row : rows) {
            summary_material += row.first;
            summary_material.push_back('|');
            summary_material += std::to_string(frontier_score_bits(row.second));
            summary_material.push_back(';');
        }
        const uint64_t summary = fnv1a64(summary_material);
        for (size_t index = 0; index < rows.size(); ++index) {
            *frontier_audit_output << event << '\t' << reads_arrived_ << '\t' << label
                << '\t' << index + 1U << '\t' << rows[index].first << '\t'
                << std::setprecision(17) << rows[index].second << '\t'
                << std::hex << std::setw(16) << std::setfill('0')
                << frontier_score_bits(rows[index].second) << std::dec << std::setfill(' ')
                << '\t' << rows[index].first << '\t' << std::hex << std::setw(16)
                << std::setfill('0') << summary << std::dec << std::setfill(' ') << '\n';
        }
    }

    void audit_candidate_frontier(const char* label, const std::vector<Candidate>& values) const {
        std::vector<std::pair<std::string, double>> rows;
        rows.reserve(values.size());
        for (const Candidate& value : values) rows.push_back({value.sequence, value.score});
        audit_frontier_rows(label, rows);
    }

    void audit_sparse_frontier(const char* label) const {
        std::vector<std::pair<std::string, double>> rows;
        rows.reserve(sparse_pool_.size());
        for (const SparseCandidateState& value : sparse_pool_) {
            rows.push_back({value.sequence, value.score});
        }
        audit_frontier_rows(label, rows);
    }

    void audit_path_frontier(const char* label, const std::vector<PathCandidate>& values) const {
        std::vector<std::pair<std::string, double>> rows;
        rows.reserve(values.size());
        for (const PathCandidate& value : values) {
            rows.push_back({materialize_path(path_arena_, value.path, config_.k), value.score});
        }
        audit_frontier_rows(label, rows);
    }

    void audit_beam_frontier(const char* label, const std::vector<BeamState>& values) const {
        std::vector<std::pair<std::string, double>> rows;
        rows.reserve(values.size());
        for (const BeamState& value : values) {
            rows.push_back({materialize_path(path_arena_, value.path, config_.k), value.cumulative_score});
        }
        audit_frontier_rows(label, rows);
    }
#endif

    void add_length(int length) {
        if (lower_lengths_.empty() || length <= lower_lengths_.top()) {
            lower_lengths_.push(length);
        } else {
            upper_lengths_.push(length);
        }
        if (lower_lengths_.size() > upper_lengths_.size() + 1U) {
            upper_lengths_.push(lower_lengths_.top());
            lower_lengths_.pop();
        } else if (upper_lengths_.size() > lower_lengths_.size()) {
            lower_lengths_.push(upper_lengths_.top());
            upper_lengths_.pop();
        }
    }

    double current_median() const {
        if (lower_lengths_.empty()) return 0.0;
        if (lower_lengths_.size() == upper_lengths_.size()) {
            return 0.5 * static_cast<double>(lower_lengths_.top() + upper_lengths_.top());
        }
        return static_cast<double>(lower_lengths_.top());
    }

    void begin_dirty_epoch() {
        dirty_contexts_.clear();
        ++current_update_epoch_;
        if (current_update_epoch_ == 0U) {
            std::fill(dirty_context_epoch_.begin(), dirty_context_epoch_.end(), 0U);
            std::fill(
                transition_score_cache_epoch_.begin(), transition_score_cache_epoch_.end(),
                std::numeric_limits<uint32_t>::max());
            std::fill(
                transition_score_cache_.begin(), transition_score_cache_.end(),
                -std::numeric_limits<double>::infinity());
            current_update_epoch_ = 1U;
        }
    }

    void mark_dirty_context(uint32_t context) {
        if (dirty_context_epoch_[context] == current_update_epoch_) return;
        dirty_context_epoch_[context] = current_update_epoch_;
        dirty_contexts_.push_back(context);
    }

    std::string decode_kmer(uint32_t code) const {
        std::string result(static_cast<size_t>(config_.k), 'A');
        for (int index = config_.k - 1; index >= 0; --index) {
            result[static_cast<size_t>(index)] = code_base(static_cast<int>(code & 3U));
            code >>= 2U;
        }
        return result;
    }

    uint32_t start_path(uint32_t code) {
        const auto found = start_paths_.find(code);
        if (found != start_paths_.end()) {
            ++path_cache_hits_;
            return found->second;
        }
        path_arena_.emplace_back();
        const uint32_t path = static_cast<uint32_t>(path_arena_.size() - 1);
        PathNode& node = path_arena_.back();
        node.parent = NO_PATH;
        node.start_code = code;
        node.context = code;
        node.length = static_cast<uint32_t>(config_.k);
        node.transitions = 0;
        node.appended_base = 0;
        start_paths_.emplace(code, path);
        ++path_nodes_created_;
        return path;
    }

    uint32_t extend_path(uint32_t parent, int base, uint32_t context) {
        const uint32_t cached = path_arena_[parent].children[static_cast<size_t>(base)];
        if (cached != NO_PATH) {
            ++path_cache_hits_;
            return cached;
        }
        path_arena_.emplace_back();
        const uint32_t path = static_cast<uint32_t>(path_arena_.size() - 1);
        path_arena_[parent].children[static_cast<size_t>(base)] = path;
        PathNode& node = path_arena_.back();
        node.parent = parent;
        node.start_code = path_arena_[parent].start_code;
        node.context = context;
        node.length = path_arena_[parent].length + 1;
        node.transitions = path_arena_[parent].transitions + 1;
        node.appended_base = static_cast<uint8_t>(base);
        ++path_nodes_created_;
        return path;
    }

    void update_graph(const std::string& sequence) {
        begin_dirty_epoch();
        uint32_t code = 0;
        int valid_run = 0;
        bool have_first = false;
        bool have_previous = false;
        uint32_t first_kmer = 0;
        uint32_t last_kmer = 0;
        uint32_t previous_kmer = 0;
        for (char ch : sequence) {
            int base = base_code(ch);
            if (base < 0) {
                code = 0;
                valid_run = 0;
                have_previous = false;
                continue;
            }
            code = ((code << 2U) | static_cast<uint32_t>(base)) & mask_;
            ++valid_run;
            if (valid_run >= config_.k) {
                node_counts_.increment(code);
                if (!have_first) {
                    first_kmer = code;
                    have_first = true;
                }
                last_kmer = code;
                if (have_previous) {
                    const uint32_t edge = (previous_kmer << 2U) | static_cast<uint32_t>(base);
                    edge_counts_.increment(edge);
                    outgoing_totals_.increment(previous_kmer);
                    mark_dirty_context(previous_kmer);
                }
                previous_kmer = code;
                have_previous = true;
            }
        }
        if (have_first) {
            start_counts_.increment(first_kmer);
            end_counts_.increment(last_kmer);
            if (total_start_observations_ == std::numeric_limits<uint64_t>::max()) {
                throw std::runtime_error("start observation counter overflow");
            }
            ++total_start_observations_;
            mark_dirty_context(first_kmer);
            mark_dirty_context(last_kmer);
        }
    }

    bool make_sparse_candidate(const std::string& sequence, SparseCandidateState& state) const {
        if (sequence.size() < static_cast<size_t>(config_.k)) return false;
        state = {};
        state.sequence = sequence;
        uint32_t code = 0;
        int valid_run = 0;
        bool have_previous = false;
        uint32_t previous = 0;
        uint32_t transition_position = 0;
        for (char ch : sequence) {
            const int base = base_code(ch);
            if (base < 0) return false;
            code = ((code << 2U) | static_cast<uint32_t>(base)) & mask_;
            ++valid_run;
            if (valid_run < config_.k) continue;
            if (!have_previous) {
                state.start_code = code;
                have_previous = true;
            } else {
                const uint32_t context = previous;
                const uint32_t edge = (context << 2U) | static_cast<uint32_t>(base);
                auto& counts = state.edge_occurrences[context];
                if (counts[static_cast<size_t>(base)] == std::numeric_limits<uint32_t>::max()) {
                    throw std::runtime_error("sparse candidate edge occurrence overflow");
                }
                ++counts[static_cast<size_t>(base)];
                state.context_positions[context].push_back(transition_position);
                const double score = transition_score(context, base);
                if (!std::isfinite(score)) return false;
                state.edge_score_cache[edge] = score;
                state.transition_sum += score;
                ++state.transitions;
                ++transition_position;
            }
            previous = code;
            state.end_code = code;
        }
        return have_previous;
    }

    double sparse_candidate_score(const SparseCandidateState& state) const {
        const auto& starts = start_counts_.empty() ? node_counts_ : start_counts_;
        const uint32_t start_support = starts.get(state.start_code);
        const double denominator = static_cast<double>(std::max<uint64_t>(1, total_start_observations_))
            + static_cast<double>(std::max<size_t>(1, starts.size()));
        const double start_probability = (static_cast<double>(start_support) + 1.0) / denominator;
        const double start_score = std::log(start_probability)
            + 0.15 * std::log1p(static_cast<double>(start_support));
        const double normalized = (start_score + state.transition_sum)
            / static_cast<double>(std::max(1, state.transitions + 1));
        return normalized
            + length_score(static_cast<int>(state.sequence.size()), current_median())
            + 0.05 * std::log1p(static_cast<double>(end_counts_.get(state.end_code)));
    }

    int refresh_sparse_candidate(SparseCandidateState& state) {
        int dependencies = 0;
        for (uint32_t context : dirty_contexts_) {
            const auto found = state.edge_occurrences.find(context);
            if (found == state.edge_occurrences.end()) continue;
            for (int base = 0; base < 4; ++base) {
                const uint32_t occurrences = found->second[static_cast<size_t>(base)];
                if (occurrences == 0) continue;
                const uint32_t edge = (context << 2U) | static_cast<uint32_t>(base);
                const double previous = state.edge_score_cache.at(edge);
                const double current = transition_score(context, base);
                if (!std::isfinite(previous) || !std::isfinite(current)) {
                    throw std::runtime_error("non-finite sparse dependency score");
                }
                state.transition_sum += static_cast<double>(occurrences) * (current - previous);
                state.edge_score_cache[edge] = current;
                dependencies += occurrences;
            }
        }
        state.score = sparse_candidate_score(state);
        return dependencies;
    }

    ConsensusResult recompute_sparse(const std::string& arriving_sequence) {
        ConsensusResult result;
        result.median_read_length = current_median();
        const auto bounds = length_bounds(result.median_read_length);
        result.minimum_length = bounds.first;
        result.maximum_length = bounds.second;
        result.dirty_context_count = static_cast<int>(dirty_contexts_.size());

        std::vector<SparseCandidateState> repaired;
        repaired.reserve(sparse_pool_.size());
        for (SparseCandidateState& candidate : sparse_pool_) {
            result.sparse_dependencies_updated += refresh_sparse_candidate(candidate);
            ++result.sparse_candidates_refreshed;

            double best_gain = 0.0;
            size_t best_index = 0;
            int best_base = -1;
            for (uint32_t context : dirty_contexts_) {
                const auto positions = candidate.context_positions.find(context);
                if (positions == candidate.context_positions.end()) continue;
                for (uint16_t transition_position : positions->second) {
                    const size_t sequence_index = static_cast<size_t>(config_.k) + transition_position;
                    if (sequence_index >= candidate.sequence.size()) continue;
                    const int old_base = base_code(candidate.sequence[sequence_index]);
                    const double old_score = transition_score(context, old_base);
                    for (int base = 0; base < 4; ++base) {
                        if (base == old_base) continue;
                        const double next_score = transition_score(context, base);
                        if (!std::isfinite(next_score)) continue;
                        const double gain = next_score - old_score;
                        if (gain > best_gain + 1e-12 ||
                            (std::fabs(gain - best_gain) <= 1e-12 && gain > 0.0 &&
                             (best_base < 0 || base < best_base))) {
                            best_gain = gain;
                            best_index = sequence_index;
                            best_base = base;
                        }
                    }
                }
            }
            if (best_base >= 0) {
                std::string repaired_sequence = candidate.sequence;
                repaired_sequence[best_index] = code_base(best_base);
                SparseCandidateState repair;
                if (make_sparse_candidate(repaired_sequence, repair)) {
                    repair.score = sparse_candidate_score(repair);
                    repaired.push_back(std::move(repair));
                    ++result.sparse_repairs_generated;
                }
            }
        }

        const bool already_present = std::any_of(
            sparse_pool_.begin(), sparse_pool_.end(), [&](const SparseCandidateState& candidate) {
                return candidate.sequence == arriving_sequence;
            });
        if (!already_present) {
            SparseCandidateState arriving;
            if (make_sparse_candidate(arriving_sequence, arriving)) {
                arriving.score = sparse_candidate_score(arriving);
                sparse_pool_.push_back(std::move(arriving));
            }
        }
        sparse_pool_.insert(
            sparse_pool_.end(),
            std::make_move_iterator(repaired.begin()),
            std::make_move_iterator(repaired.end()));

        std::sort(sparse_pool_.begin(), sparse_pool_.end(), [](const auto& lhs, const auto& rhs) {
            const int score_order = compare_finite_scores_descending(lhs.score, rhs.score);
            if (score_order != 0) return score_order < 0;
            return lhs.sequence < rhs.sequence;
        });
        sparse_pool_.erase(
            std::unique(sparse_pool_.begin(), sparse_pool_.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.sequence == rhs.sequence;
            }),
            sparse_pool_.end());
        if (sparse_pool_.size() > static_cast<size_t>(config_.beam_width)) {
            sparse_pool_.resize(static_cast<size_t>(config_.beam_width));
        }
#ifdef ACOR_FRONTIER_AUDIT
        audit_sparse_frontier("sparse_pool_prune");
#endif
        for (const SparseCandidateState& candidate : sparse_pool_) {
            result.ranked.push_back({candidate.sequence, candidate.score});
        }
        if (result.ranked.empty()) {
            const std::string fallback = fallback_read();
            if (!fallback.empty()) result.ranked.push_back({fallback, 0.0});
        }
        result.selection_mode = "sparse_dirty_candidate_pool";
        if (!result.ranked.empty()) {
            result.selected_length_delta =
                static_cast<int>(result.ranked.front().sequence.size()) - design_length_;
        }
        sparse_candidates_refreshed_total_ += result.sparse_candidates_refreshed;
        sparse_dependencies_updated_total_ += result.sparse_dependencies_updated;
        sparse_repairs_generated_total_ += result.sparse_repairs_generated;
        return result;
    }

    std::string fallback_read() const {
        std::pair<std::string, uint32_t> best{"", 0};
        for (const auto& item : read_counts_) {
            if (item.second > best.second || (item.second == best.second && item.first < best.first)) best = item;
        }
        return best.first;
    }

    std::string fallback_read_exact_length() const {
        std::pair<std::string, uint32_t> best{"", 0};
        for (const auto& item : read_counts_) {
            if (static_cast<int>(item.first.size()) != design_length_) continue;
            if (item.second > best.second || (item.second == best.second && item.first < best.first)) best = item;
        }
        return best.first;
    }

    std::pair<int, int> length_bounds(double med) const {
        if (config_.length_prior == "median_band") {
            const int tolerance = std::max(3, static_cast<int>(std::lround(0.15 * med)));
            return {std::max(config_.k, static_cast<int>(std::floor(med)) - tolerance),
                    std::max(config_.k, static_cast<int>(std::ceil(med)) + tolerance)};
        }
        const double count = static_cast<double>(std::max<size_t>(1, reads_arrived_));
        const double mu = log_length_sum_ / count;
        const double variance = std::max(0.0, log_length_sq_sum_ / count - mu * mu);
        const double sigma = std::max(0.08, std::sqrt(variance));
        int minimum = static_cast<int>(std::floor(std::exp(mu - 2.0 * sigma)));
        int maximum = static_cast<int>(std::ceil(std::exp(mu + 2.0 * sigma)));
        minimum = std::max(config_.k, std::max(minimum, static_cast<int>(std::floor(0.60 * med))));
        maximum = std::max(minimum, std::min(maximum, static_cast<int>(std::ceil(1.50 * med))));
        return {minimum, maximum};
    }

    double length_score(int length, double med) const {
        if (config_.length_prior == "median_band") {
            const double tolerance = static_cast<double>(std::max(3, static_cast<int>(std::lround(0.15 * med))));
            return -0.75 * std::fabs(static_cast<double>(length) - med) / tolerance;
        }
        const double count = static_cast<double>(std::max<size_t>(1, reads_arrived_));
        const double mu = log_length_sum_ / count;
        const double variance = std::max(0.0, log_length_sq_sum_ / count - mu * mu);
        const double sigma = std::max(0.08, std::sqrt(variance));
        const double z = (std::log(static_cast<double>(std::max(1, length))) - mu) / sigma;
        return -0.375 * z * z;
    }

    void retain_top_candidates(std::vector<Candidate>& candidates) const {
        std::sort(candidates.begin(), candidates.end(), candidate_better);
        candidates.erase(std::unique(candidates.begin(), candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
            return lhs.sequence == rhs.sequence;
        }), candidates.end());
        if (candidates.size() > static_cast<size_t>(config_.beam_width)) {
            candidates.resize(static_cast<size_t>(config_.beam_width));
        }
#ifdef ACOR_FRONTIER_AUDIT
        audit_candidate_frontier("candidate_prune", candidates);
#endif
    }

    void retain_top_path_candidates(std::vector<PathCandidate>& candidates) const {
        std::sort(candidates.begin(), candidates.end(), [this](const PathCandidate& lhs, const PathCandidate& rhs) {
            const int score_order = compare_finite_scores_descending(lhs.score, rhs.score);
            if (score_order != 0) return score_order < 0;
            if (path_arena_[lhs.path].length == path_arena_[rhs.path].length) {
                return path_lex_less_equal_length(path_arena_, lhs.path, rhs.path);
            }
            return materialize_path(path_arena_, lhs.path, config_.k)
                < materialize_path(path_arena_, rhs.path, config_.k);
        });
        candidates.erase(std::unique(candidates.begin(), candidates.end(), [](const PathCandidate& lhs, const PathCandidate& rhs) {
            return lhs.path == rhs.path;
        }), candidates.end());
        if (candidates.size() > static_cast<size_t>(config_.beam_width)) {
            candidates.resize(static_cast<size_t>(config_.beam_width));
        }
#ifdef ACOR_FRONTIER_AUDIT
        audit_path_frontier("path_candidate_prune", candidates);
#endif
    }

    std::vector<Candidate> materialize(const std::vector<PathCandidate>& candidates) const {
        std::vector<Candidate> result;
        result.reserve(candidates.size());
        for (const PathCandidate& candidate : candidates) {
            result.push_back({materialize_path(path_arena_, candidate.path, config_.k), candidate.score});
        }
        return result;
    }

    double transition_score(uint32_t context, int base) const {
        if (transition_score_cache_epoch_[context] != dirty_context_epoch_[context]) {
            const uint32_t total = outgoing_totals_.get(context);
            for (int candidate_base = 0; candidate_base < 4; ++candidate_base) {
                const uint32_t candidate_edge =
                    (context << 2U) | static_cast<uint32_t>(candidate_base);
                const uint32_t support = edge_counts_.get(candidate_edge);
                transition_score_cache_[candidate_edge] = support == 0U
                    ? -std::numeric_limits<double>::infinity()
                    : std::log(
                        (static_cast<double>(support) + 1.0) /
                        (static_cast<double>(total) + 4.0))
                        + 0.15 * std::log1p(static_cast<double>(support));
            }
            transition_score_cache_epoch_[context] = dirty_context_epoch_[context];
        }
        const uint32_t edge = (context << 2U) | static_cast<uint32_t>(base);
        return transition_score_cache_[edge];
    }

    BandPlan make_band_plan(int target_length) const {
        const int transitions = std::max(0, target_length - config_.k);
        if (optimization_mode_ != OptimizationMode::Banded || band_width_ >= transitions) {
            return {0, transitions, transitions / 2};
        }
        const auto& starts = start_counts_.empty() ? node_counts_ : start_counts_;
        uint32_t context = 0;
        uint32_t best_start_support = 0;
        bool have_start = false;
        for (uint32_t code : starts.touched()) {
            const uint32_t support = starts.get(code);
            if (!have_start || support > best_start_support ||
                (support == best_start_support && code < context)) {
                context = code;
                best_start_support = support;
                have_start = true;
            }
        }
        int weakest = 0;
        double weakest_margin = std::numeric_limits<double>::infinity();
        for (int transition = 0; have_start && transition < transitions; ++transition) {
            double best = -std::numeric_limits<double>::infinity();
            double second = -std::numeric_limits<double>::infinity();
            int best_base = 0;
            for (int base = 0; base < 4; ++base) {
                const double score = transition_score(context, base);
                if (score > best + 1e-12 || (std::fabs(score - best) <= 1e-12 && base < best_base)) {
                    second = best;
                    best = score;
                    best_base = base;
                } else if (score > second) {
                    second = score;
                }
            }
            if (!std::isfinite(best)) break;
            const double margin = std::isfinite(second)
                ? best - second : std::numeric_limits<double>::infinity();
            if (margin < weakest_margin) {
                weakest_margin = margin;
                weakest = transition;
            }
            context = ((context << 2U) | static_cast<uint32_t>(best_base)) & mask_;
        }
        int start = std::max(0, weakest - band_width_ / 2);
        int end = std::min(transitions, start + band_width_);
        start = std::max(0, end - band_width_);
        return {start, end, weakest};
    }

    size_t beam_limit(int transition, const BandPlan& plan) const {
        if (optimization_mode_ != OptimizationMode::Banded) {
            return static_cast<size_t>(config_.beam_width);
        }
        return transition >= plan.start_transition && transition < plan.end_transition
            ? static_cast<size_t>(config_.beam_width) : 1U;
    }

    ConsensusResult recompute() {
        ConsensusResult result;
        result.median_read_length = current_median();
        auto bounds = length_bounds(result.median_read_length);
        if (constraint_mode_ == ConstraintMode::HardEndpoint) {
            bounds = {design_length_, design_length_};
        }
        result.minimum_length = bounds.first;
        result.maximum_length = bounds.second;
        const BandPlan band = make_band_plan(result.maximum_length);
        result.band_start_transition = band.start_transition;
        result.band_end_transition = band.end_transition;
        result.band_center_transition = band.center_transition;

        const auto& starts = start_counts_.empty() ? node_counts_ : start_counts_;
        std::vector<std::pair<uint32_t, uint32_t>> ranked_starts;
        ranked_starts.reserve(starts.size());
        for (uint32_t code : starts.touched()) ranked_starts.push_back({code, starts.get(code)});
        std::sort(ranked_starts.begin(), ranked_starts.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.second != rhs.second) return lhs.second > rhs.second;
            return lhs.first < rhs.first;
        });
        const size_t start_limit = beam_limit(0, band);
        if (ranked_starts.size() > start_limit) {
            ranked_starts.resize(start_limit);
        }
        uint64_t total_starts = 0;
        for (uint32_t code : starts.touched()) total_starts += starts.get(code);

        std::vector<BeamState> beam;
        for (const auto& item : ranked_starts) {
            const double probability = (static_cast<double>(item.second) + 1.0) /
                (static_cast<double>(total_starts) + static_cast<double>(std::max<size_t>(1, starts.size())));
            beam.push_back({start_path(item.first),
                            std::log(probability) + 0.15 * std::log1p(static_cast<double>(item.second))});
        }

        std::vector<PathCandidate> complete;
        while (!beam.empty()) {
            const int current_length = static_cast<int>(path_arena_[beam.front().path].length);
            if (current_length >= result.minimum_length) {
                for (const BeamState& state : beam) {
                    const PathNode& path = path_arena_[state.path];
                    const uint32_t end_support = end_counts_.get(path.context);
                    const double normalized = state.cumulative_score /
                        static_cast<double>(std::max(1, path.transitions + 1));
                    complete.push_back({state.path, normalized + length_score(current_length, result.median_read_length)
                                        + 0.05 * std::log1p(static_cast<double>(end_support))});
                }
                retain_top_path_candidates(complete);
            }
            if (current_length >= result.maximum_length) break;

            std::vector<BeamState> next;
            next.reserve(beam.size() * 4U);
            for (const BeamState& state : beam) {
                const uint32_t context = path_arena_[state.path].context;
                const uint32_t total = outgoing_totals_.get(context);
                for (int base = 0; base < 4; ++base) {
                    const uint32_t edge = (context << 2U) | static_cast<uint32_t>(base);
                    const uint32_t support = edge_counts_.get(edge);
                    if (support == 0U) continue;
                    const double probability = (static_cast<double>(support) + 1.0) /
                        (static_cast<double>(total) + 4.0);
                    BeamState extension{
                        extend_path(state.path, base, edge & mask_),
                        state.cumulative_score + (std::log(probability) +
                            0.15 * std::log1p(static_cast<double>(support))),
                    };
                    next.push_back(std::move(extension));
                }
            }
            std::sort(next.begin(), next.end(), [this](const BeamState& lhs, const BeamState& rhs) {
                const int score_order = compare_finite_scores_descending(
                    lhs.cumulative_score, rhs.cumulative_score);
                if (score_order != 0) return score_order < 0;
                return path_lex_less_equal_length(path_arena_, lhs.path, rhs.path);
            });
            next.erase(std::unique(next.begin(), next.end(), [](const BeamState& lhs, const BeamState& rhs) {
                return lhs.path == rhs.path;
            }), next.end());
            const int next_transition = current_length - config_.k;
            const size_t limit = beam_limit(next_transition, band);
            if (next.size() > limit) {
                next.resize(limit);
            }
#ifdef ACOR_FRONTIER_AUDIT
            audit_beam_frontier("recompute_beam_prune", next);
#endif
            beam = std::move(next);
        }

        retain_top_path_candidates(complete);
        const bool selected_beam_candidate = !complete.empty();
        result.exact_m_candidate_count = constraint_mode_ == ConstraintMode::HardEndpoint
            ? static_cast<int>(complete.size()) : 0;
        result.ranked = materialize(complete);
        if (result.ranked.empty()) {
            const std::string fallback = constraint_mode_ == ConstraintMode::HardEndpoint
                ? fallback_read_exact_length() : fallback_read();
            if (!fallback.empty()) result.ranked.push_back({fallback, 0.0});
        }
        result.selection_mode = constraint_mode_ == ConstraintMode::HardEndpoint
            ? (selected_beam_candidate ? "exact_m_beam" : "exact_m_read_fallback")
            : "arm0_consensus";
        if (!result.ranked.empty()) {
            result.selected_length_delta = static_cast<int>(result.ranked.front().sequence.size()) - design_length_;
        }
        return result;
    }

    ConsensusResult recompute_soft() {
        ConsensusResult result;
        result.median_read_length = current_median();
        const auto arm0_bounds = length_bounds(result.median_read_length);
        const int soft_min = std::max(config_.k, design_length_ - 5);
        const int soft_max = std::max(soft_min, design_length_ + 5);
        result.minimum_length = soft_min;
        result.maximum_length = soft_max;
        const int search_max = std::max(arm0_bounds.second, soft_max);
        const BandPlan band = make_band_plan(design_length_);
        result.band_start_transition = band.start_transition;
        result.band_end_transition = band.end_transition;
        result.band_center_transition = band.center_transition;

        const auto& starts = start_counts_.empty() ? node_counts_ : start_counts_;
        std::vector<std::pair<uint32_t, uint32_t>> ranked_starts;
        ranked_starts.reserve(starts.size());
        for (uint32_t code : starts.touched()) ranked_starts.push_back({code, starts.get(code)});
        std::sort(ranked_starts.begin(), ranked_starts.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.second != rhs.second) return lhs.second > rhs.second;
            return lhs.first < rhs.first;
        });
        const size_t start_limit = beam_limit(0, band);
        if (ranked_starts.size() > start_limit) {
            ranked_starts.resize(start_limit);
        }
        uint64_t total_starts = 0;
        for (uint32_t code : starts.touched()) total_starts += starts.get(code);

        std::vector<BeamState> beam;
        for (const auto& item : ranked_starts) {
            const double probability = (static_cast<double>(item.second) + 1.0) /
                (static_cast<double>(total_starts) + static_cast<double>(std::max<size_t>(1, starts.size())));
            beam.push_back({start_path(item.first),
                            std::log(probability) + 0.15 * std::log1p(static_cast<double>(item.second))});
        }

        std::vector<PathCandidate> arm0_complete;
        std::vector<std::vector<PathCandidate>> soft_by_length(static_cast<size_t>(soft_max - soft_min + 1));
        while (!beam.empty()) {
            const int current_length = static_cast<int>(path_arena_[beam.front().path].length);
            if (current_length >= arm0_bounds.first && current_length <= arm0_bounds.second) {
                for (const BeamState& state : beam) {
                    const PathNode& path = path_arena_[state.path];
                    const uint32_t end_support = end_counts_.get(path.context);
                    const double normalized = state.cumulative_score /
                        static_cast<double>(std::max(1, path.transitions + 1));
                    arm0_complete.push_back({state.path,
                        normalized + length_score(current_length, result.median_read_length)
                        + 0.05 * std::log1p(static_cast<double>(end_support))});
                }
                retain_top_path_candidates(arm0_complete);
            }
            if (current_length >= soft_min && current_length <= soft_max) {
                auto& at_length = soft_by_length[static_cast<size_t>(current_length - soft_min)];
                for (const BeamState& state : beam) {
                    const PathNode& path = path_arena_[state.path];
                    const uint32_t end_support = end_counts_.get(path.context);
                    const double normalized = state.cumulative_score /
                        static_cast<double>(std::max(1, path.transitions + 1));
                    at_length.push_back({state.path,
                        normalized + length_score(current_length, result.median_read_length)
                        + 0.05 * std::log1p(static_cast<double>(end_support))});
                }
                retain_top_path_candidates(at_length);
            }
            if (current_length >= search_max) break;

            std::vector<BeamState> next;
            next.reserve(beam.size() * 4U);
            for (const BeamState& state : beam) {
                const uint32_t context = path_arena_[state.path].context;
                const uint32_t total = outgoing_totals_.get(context);
                for (int base = 0; base < 4; ++base) {
                    const uint32_t edge = (context << 2U) | static_cast<uint32_t>(base);
                    const uint32_t support = edge_counts_.get(edge);
                    if (support == 0U) continue;
                    const double probability = (static_cast<double>(support) + 1.0) /
                        (static_cast<double>(total) + 4.0);
                    BeamState extension{
                        extend_path(state.path, base, edge & mask_),
                        state.cumulative_score + (std::log(probability) +
                            0.15 * std::log1p(static_cast<double>(support))),
                    };
                    next.push_back(std::move(extension));
                }
            }
            std::sort(next.begin(), next.end(), [this](const BeamState& lhs, const BeamState& rhs) {
                const int score_order = compare_finite_scores_descending(
                    lhs.cumulative_score, rhs.cumulative_score);
                if (score_order != 0) return score_order < 0;
                return path_lex_less_equal_length(path_arena_, lhs.path, rhs.path);
            });
            next.erase(std::unique(next.begin(), next.end(), [](const BeamState& lhs, const BeamState& rhs) {
                return lhs.path == rhs.path;
            }), next.end());
            const int next_transition = current_length - config_.k;
            const size_t limit = beam_limit(next_transition, band);
            if (next.size() > limit) {
                next.resize(limit);
            }
#ifdef ACOR_FRONTIER_AUDIT
            audit_beam_frontier("soft_beam_prune", next);
#endif
            beam = std::move(next);
        }

        for (const auto& candidates : soft_by_length) {
            result.soft_window_candidate_count += static_cast<int>(candidates.size());
        }
        std::vector<PathCandidate> exact_m;
        if (design_length_ >= soft_min && design_length_ <= soft_max) {
            exact_m = soft_by_length[static_cast<size_t>(design_length_ - soft_min)];
        }
        result.exact_m_candidate_count = static_cast<int>(exact_m.size());
        if (exact_m.empty()) {
            const std::string exact_read = fallback_read_exact_length();
            if (!exact_read.empty()) {
                result.ranked.push_back({exact_read, 0.0});
                result.selection_mode = "exact_m_read_fallback";
                result.selected_length_delta = 0;
                return result;
            }
        }
        if (!exact_m.empty()) {
            retain_top_path_candidates(exact_m);
            result.ranked = materialize(exact_m);
            result.selection_mode = "exact_m_beam";
            result.selected_length_delta = 0;
            return result;
        }

        for (int distance = 1; distance <= 5; ++distance) {
            std::vector<PathCandidate> nearest;
            for (int target : {design_length_ - distance, design_length_ + distance}) {
                if (target < soft_min || target > soft_max) continue;
                const auto& candidates = soft_by_length[static_cast<size_t>(target - soft_min)];
                nearest.insert(nearest.end(), candidates.begin(), candidates.end());
            }
            if (!nearest.empty()) {
                retain_top_path_candidates(nearest);
                result.ranked = materialize(nearest);
                result.selection_mode = "soft_nearest_beam";
                result.selected_length_delta =
                    static_cast<int>(result.ranked.front().sequence.size()) - design_length_;
                return result;
            }
        }

        retain_top_path_candidates(arm0_complete);
        result.ranked = materialize(arm0_complete);
        if (result.ranked.empty()) {
            const std::string fallback = fallback_read();
            if (!fallback.empty()) result.ranked.push_back({fallback, 0.0});
        }
        result.selection_mode = "arm0_fallback";
        result.fallback_to_arm0 = true;
        if (!result.ranked.empty()) {
            result.selected_length_delta = static_cast<int>(result.ranked.front().sequence.size()) - design_length_;
        }
        return result;
    }
};

struct AcorKernel::Impl {
    Config config = parse_config("kmc_k9_b16_lognormal");
    IncrementalKmcBeam model{
        config, ConstraintMode::SoftLength, OptimizationMode::FinalOnlyExact, 0};
};

AcorKernel::AcorKernel() : impl_(std::make_unique<Impl>()) {}
AcorKernel::~AcorKernel() = default;
AcorKernel::AcorKernel(AcorKernel&&) noexcept = default;
AcorKernel& AcorKernel::operator=(AcorKernel&&) noexcept = default;

void AcorKernel::reset(int design_length) {
    impl_->model.reset(design_length);
}

AcorIncrementalResult AcorKernel::add_read(const std::string& sequence, bool is_last) {
    const ConsensusResult observed = impl_->model.add_read(sequence, is_last);
    AcorIncrementalResult output;
    output.final_decode_seconds = observed.final_decode_seconds;
    output.ranked.reserve(observed.ranked.size());
    for (const Candidate& candidate : observed.ranked) {
        output.ranked.push_back({candidate.sequence, candidate.score});
    }
    return output;
}

double AcorKernel::graph_update_seconds() const {
    return impl_->model.graph_update_seconds();
}

double AcorKernel::beam_decode_seconds() const {
    return impl_->model.beam_decode_seconds();
}

std::string acor_canonical(const std::string& sequence) {
    return canonical(sequence);
}

void acor_edit_distance_self_test() {
    edit_distance_self_test();
}
