#include "acor_kernel.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <pthread.h>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <sys/resource.h>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;

struct ClusterTask {
    std::string cluster_id;
    std::vector<std::string> reads;
    int design_length = 0;
};

struct ClusterResult {
    std::string cluster_id;
    std::string candidate;
    size_t stop_index = 0;
    size_t total_reads = 0;
    size_t legacy_prefix_alias = 0;
    bool stopped_early = false;
    bool completed_empty = false;
    size_t worker_id = 0;
    double model_reset_seconds = 0.0;
    double canonical_seconds = 0.0;
    double add_read_seconds = 0.0;
    double soft_decode_seconds = 0.0;
    double graph_update_seconds = 0.0;
    double beam_decode_seconds = 0.0;
    double evidence_ed_seconds = 0.0;
    double cluster_seconds = 0.0;
    double ledger_accumulation_seconds = 0.0;
    double ledger_merge_seconds = 0.0;
    double ledger_final_decode_seconds = 0.0;
};

struct WorkerStats {
    size_t clusters = 0;
    size_t reads_consumed = 0;
    double cluster_seconds_sum = 0.0;
};

struct WorkerContext {
    AcorKernel model;
    bool hot_cache_enabled = false;
};

struct LedgerBlockTask {
    size_t cluster_index = 0;
    size_t block_index = 0;
    size_t begin = 0;
    size_t end = 0;
};

struct LedgerBlockResult {
    size_t cluster_index = 0;
    size_t block_index = 0;
    size_t worker_id = 0;
    size_t reads = 0;
    double canonical_seconds = 0.0;
    double evidence_seconds = 0.0;
    double block_seconds = 0.0;
};

struct SparseEvidenceBlock {
    std::unordered_map<std::uint64_t, std::uint32_t> counters;
    std::unordered_map<std::string, std::uint32_t> read_counts;
    std::vector<int> lengths;
    std::size_t reads_arrived = 0;
    std::size_t exact_m_read_count = 0;
    double log_length_sum = 0.0;
    double log_length_sq_sum = 0.0;
    std::uint64_t total_start_observations = 0;

    static std::uint64_t key(std::uint8_t kind, std::uint32_t index) {
        return (static_cast<std::uint64_t>(kind) << 32U) | index;
    }

    void add_counter(std::uint8_t kind, std::uint32_t index) {
        const std::uint64_t encoded = key(kind, index);
        std::uint32_t& target = counters[encoded];
        if (target == std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("sparse counter overflow");
        }
        ++target;
    }

    void add_read(const std::string& sequence, int design_length, int k) {
        if (reads_arrived == std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("sparse read count overflow");
        }
        ++reads_arrived;
        std::uint32_t& multiplicity = read_counts[sequence];
        if (multiplicity == std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("sparse read multiplicity overflow");
        }
        ++multiplicity;
        lengths.push_back(static_cast<int>(sequence.size()));
        const double log_length =
            std::log(static_cast<double>(std::max<std::size_t>(1, sequence.size())));
        log_length_sum += log_length;
        log_length_sq_sum += log_length * log_length;
        if (static_cast<int>(sequence.size()) == design_length) ++exact_m_read_count;

        const std::uint32_t mask = (static_cast<std::uint32_t>(1U) << (2 * k)) - 1U;
        std::uint32_t code = 0U;
        std::uint32_t previous = 0U;
        std::uint32_t first = 0U;
        std::uint32_t last = 0U;
        int valid_run = 0;
        bool have_previous = false;
        bool have_first = false;
        for (char ch : sequence) {
            int base = -1;
            switch (ch) {
            case 'A': base = 0; break;
            case 'C': base = 1; break;
            case 'G': base = 2; break;
            case 'T': base = 3; break;
            default: break;
            }
            if (base < 0) {
                code = 0U;
                valid_run = 0;
                have_previous = false;
                continue;
            }
            code = ((code << 2U) | static_cast<std::uint32_t>(base)) & mask;
            ++valid_run;
            if (valid_run < k) continue;
            add_counter(0U, code);
            if (!have_first) {
                first = code;
                have_first = true;
            }
            last = code;
            if (have_previous) {
                add_counter(1U, (previous << 2U) | static_cast<std::uint32_t>(base));
                add_counter(2U, previous);
            }
            previous = code;
            have_previous = true;
        }
        if (have_first) {
            add_counter(3U, first);
            add_counter(4U, last);
            ++total_start_observations;
        }
    }
};

// Metadata-only companion for the shared-atomic experiment.  Graph counters
// are updated directly in one shared atomic ledger; read multiplicities and
// length statistics remain private to workers and are merged deterministically
// after the counter phase.
struct AtomicMetadataBlock {
    std::unordered_map<std::string, std::uint32_t> read_counts;
    std::vector<int> lengths;
    std::size_t reads_arrived = 0;
    std::size_t exact_m_read_count = 0;
    double log_length_sum = 0.0;
    double log_length_sq_sum = 0.0;
    std::uint64_t total_start_observations = 0;

    void add_read(const std::string& sequence, int design_length, int k) {
        ++reads_arrived;
        ++read_counts[sequence];
        lengths.push_back(static_cast<int>(sequence.size()));
        const double log_length = std::log(
            static_cast<double>(std::max<std::size_t>(1U, sequence.size())));
        log_length_sum += log_length;
        log_length_sq_sum += log_length * log_length;
        if (static_cast<int>(sequence.size()) == design_length) ++exact_m_read_count;

        int valid_run = 0;
        bool have_first = false;
        for (char ch : sequence) {
            const bool valid = ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T';
            if (!valid) {
                valid_run = 0;
                continue;
            }
            ++valid_run;
            if (valid_run >= k) have_first = true;
        }
        if (have_first) ++total_start_observations;
    }
};

// Sparse tile ledger for the owner-tile experiment.  The tile directory is
// fixed and pointer-sized, while tile payloads are allocated only when a
// worker actually touches that region of the k-mer counter space.  Updates
// stay local and use a compact 256-counter payload; no atomics or per-event
// raw update objects are created.
struct TileLedger {
    static constexpr std::size_t tile_size = 256U;
    static constexpr std::size_t node_slots = 1U << 18U;
    static constexpr std::size_t edge_slots = 1U << 20U;
    static constexpr std::size_t tiles_per_node = node_slots / tile_size;
    static constexpr std::size_t tiles_per_edge = edge_slots / tile_size;
    static constexpr std::size_t total_tiles =
        tiles_per_node * 4U + tiles_per_edge;

    struct Tile {
        std::array<std::uint32_t, tile_size> counts{};
        std::vector<std::uint16_t> touched;
    };

    std::vector<std::unique_ptr<Tile>> tiles;

    TileLedger() : tiles(total_tiles) {}

    static std::size_t tile_id(std::uint8_t kind, std::uint32_t index) {
        const std::size_t base =
            kind == 0U ? 0U :
            kind == 1U ? tiles_per_node :
            kind == 2U ? tiles_per_node + tiles_per_edge :
            kind == 3U ? tiles_per_node * 2U + tiles_per_edge :
                         tiles_per_node * 3U + tiles_per_edge;
        return base + index / tile_size;
    }

    void add_counter(std::uint8_t kind, std::uint32_t index) {
        const std::size_t id = tile_id(kind, index);
        if (!tiles[id]) tiles[id] = std::make_unique<Tile>();
        Tile& tile = *tiles[id];
        const std::size_t offset = index % tile_size;
        if (tile.counts[offset] == std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("owner tile counter overflow");
        }
        if (tile.counts[offset] == 0U) {
            tile.touched.push_back(static_cast<std::uint16_t>(offset));
        }
        ++tile.counts[offset];
    }

    void add_sequence(const std::string& sequence, int k) {
        std::uint32_t code = 0U;
        std::uint32_t previous = 0U;
        std::uint32_t first = 0U;
        std::uint32_t last = 0U;
        int valid_run = 0;
        bool have_previous = false;
        bool have_first = false;
        const std::uint32_t mask = (static_cast<std::uint32_t>(1U) << (2 * k)) - 1U;
        for (char ch : sequence) {
            int base = -1;
            switch (ch) {
            case 'A': base = 0; break;
            case 'C': base = 1; break;
            case 'G': base = 2; break;
            case 'T': base = 3; break;
            default: break;
            }
            if (base < 0) {
                code = 0U;
                valid_run = 0;
                have_previous = false;
                continue;
            }
            code = ((code << 2U) | static_cast<std::uint32_t>(base)) & mask;
            ++valid_run;
            if (valid_run < k) continue;
            add_counter(0U, code);
            if (!have_first) {
                first = code;
                have_first = true;
            }
            last = code;
            if (have_previous) {
                add_counter(1U, (previous << 2U) | static_cast<std::uint32_t>(base));
                add_counter(2U, previous);
            }
            previous = code;
            have_previous = true;
        }
        if (have_first) {
            add_counter(3U, first);
            add_counter(4U, last);
        }
    }
};

// Batched owner-reduction block.  Unlike SparseEvidenceBlock, graph evidence
// is kept as a short-lived append-only update buffer.  The buffer is reduced
// by owner after each batch, so memory is bounded by batch size rather than by
// the full number of reads.
struct BatchEvidenceBlock {
    std::vector<AcorCounterDelta> deltas;
    std::unordered_map<std::string, std::uint32_t> read_counts;
    std::vector<int> lengths;
    std::size_t reads_arrived = 0;
    std::size_t exact_m_read_count = 0;
    double log_length_sum = 0.0;
    double log_length_sq_sum = 0.0;
    std::uint64_t total_start_observations = 0;
    double canonical_seconds = 0.0;
    double evidence_seconds = 0.0;

    void clear_batch() {
        deltas.clear();
        canonical_seconds = 0.0;
        evidence_seconds = 0.0;
    }

    void add_read(
        const std::string& sequence, int design_length, const AcorKernel& emitter,
        bool profile_enabled) {
        if (reads_arrived == std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("batched read count overflow");
        }
        ++reads_arrived;
        std::uint32_t& multiplicity = read_counts[sequence];
        if (multiplicity == std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("batched read multiplicity overflow");
        }
        ++multiplicity;
        lengths.push_back(static_cast<int>(sequence.size()));
        const double log_length =
            std::log(static_cast<double>(std::max<std::size_t>(1, sequence.size())));
        log_length_sum += log_length;
        log_length_sq_sum += log_length * log_length;
        if (static_cast<int>(sequence.size()) == design_length) ++exact_m_read_count;
        const std::size_t before = deltas.size();
        const auto started = std::chrono::steady_clock::now();
        emitter.append_evidence_deltas(sequence, deltas);
        if (profile_enabled) {
            evidence_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
        }
        for (std::size_t index = before; index < deltas.size(); ++index) {
            if (deltas[index].kind == 3U) ++total_start_observations;
        }
    }
};

static std::size_t sparse_owner(std::uint64_t key, std::size_t owners) {
    // SplitMix64 gives a stable, inexpensive ownership function even when
    // counter indices have regular low-bit patterns.
    key += 0x9e3779b97f4a7c15ULL;
    key = (key ^ (key >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    key = (key ^ (key >> 27U)) * 0x94d049bb133111ebULL;
    key ^= key >> 31U;
    return static_cast<std::size_t>(key % owners);
}

struct CpuLocation {
    int cpu = -1;
    int socket = -1;
    int core = -1;
};

struct OrderRecord {
    std::string cluster_id;
    std::vector<size_t> permutation;
};

class GzipReader {
public:
    explicit GzipReader(const fs::path& path) : handle_(gzopen(path.c_str(), "rb")) {
        if (handle_ == nullptr) throw std::runtime_error("cannot open gzip order pack");
        gzbuffer(handle_, 1024U * 1024U);
    }

    GzipReader(const GzipReader&) = delete;
    GzipReader& operator=(const GzipReader&) = delete;

    ~GzipReader() {
        if (handle_ != nullptr) gzclose(handle_);
    }

    gzFile get() const { return handle_; }

    void close_checked() {
        if (handle_ != nullptr && gzclose(handle_) != Z_OK) {
            handle_ = nullptr;
            throw std::runtime_error("gzip close failed");
        }
        handle_ = nullptr;
    }

private:
    gzFile handle_ = nullptr;
};

class ThreadJoiner {
public:
    explicit ThreadJoiner(std::vector<std::thread>& threads) : threads_(threads) {}
    ~ThreadJoiner() {
        for (std::thread& thread : threads_) {
            if (thread.joinable()) thread.join();
        }
    }

private:
    std::vector<std::thread>& threads_;
};

static double timeval_seconds(const timeval& value) {
    return static_cast<double>(value.tv_sec) + static_cast<double>(value.tv_usec) / 1e6;
}

// The dynamic worker pool already balances work at cluster granularity, but
// FIFO dispatch can still leave one worker with the heaviest clusters near the
// end of a run.  Sorting the independent tasks by a deterministic, inexpensive
// work estimate gives the pool a longest-processing-time-first schedule.  This
// changes only dispatch order; results remain stored by their original index,
// so reconstruction semantics and output ordering are unchanged.
static void sort_indices_by_estimated_work(
    std::vector<size_t>& indices, const std::vector<ClusterTask>& tasks) {
    std::stable_sort(
        indices.begin(), indices.end(), [&](size_t left, size_t right) {
            const size_t left_work = tasks[left].reads.size();
            const size_t right_work = tasks[right].reads.size();
            if (left_work != right_work) return left_work > right_work;
            return left < right;
        });
}

static std::string read_text_line(const fs::path& path) {
    std::ifstream handle(path);
    std::string value;
    if (!handle || !std::getline(handle, value)) {
        throw std::runtime_error("cannot read " + path.string());
    }
    return value;
}

static int read_integer_file(const fs::path& path) {
    const std::string text = read_text_line(path);
    int value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc()) {
        throw std::runtime_error("invalid integer in " + path.string());
    }
    return value;
}

static CpuLocation cpu_location(int cpu) {
    const fs::path root = fs::path("/sys/devices/system/cpu") / ("cpu" + std::to_string(cpu)) / "topology";
    return {
        cpu,
        read_integer_file(root / "physical_package_id"),
        read_integer_file(root / "core_id"),
    };
}

static void pin_current_thread(int cpu) {
    if (cpu < 0 || cpu >= CPU_SETSIZE) {
        throw std::runtime_error("cpu outside CPU_SETSIZE");
    }
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(cpu, &mask);
    const int status = pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask);
    if (status != 0) {
        throw std::runtime_error("pthread_setaffinity_np failed for cpu " + std::to_string(cpu));
    }
    cpu_set_t observed;
    CPU_ZERO(&observed);
    if (pthread_getaffinity_np(pthread_self(), sizeof(observed), &observed) != 0 ||
        CPU_COUNT(&observed) != 1 || !CPU_ISSET(cpu, &observed)) {
        throw std::runtime_error("affinity verification failed for cpu " + std::to_string(cpu));
    }
}

static std::vector<int> parse_cpu_list(const std::string& text) {
    std::vector<int> result;
    size_t begin = 0;
    while (begin <= text.size()) {
        const size_t end = text.find(',', begin);
        const std::string_view token(text.data() + begin,
                                     (end == std::string::npos ? text.size() : end) - begin);
        if (token.empty()) throw std::runtime_error("empty cpu-list token");
        int cpu = -1;
        const auto parsed = std::from_chars(token.data(), token.data() + token.size(), cpu);
        if (parsed.ec != std::errc() || parsed.ptr != token.data() + token.size() || cpu < 0) {
            throw std::runtime_error("invalid cpu-list token");
        }
        if (cpu >= CPU_SETSIZE) throw std::runtime_error("cpu-list token outside CPU_SETSIZE");
        result.push_back(cpu);
        if (end == std::string::npos) break;
        begin = end + 1U;
    }
    return result;
}

static std::string load_average() {
    try {
        return read_text_line("/proc/loadavg");
    } catch (...) {
        return "unavailable";
    }
}

static std::array<std::string_view, 7> split_input_row(const std::string& line) {
    std::array<std::string_view, 7> fields;
    size_t begin = 0;
    for (size_t index = 0; index < fields.size(); ++index) {
        const size_t end = line.find('\t', begin);
        if (index + 1U == fields.size()) {
            if (end != std::string::npos) throw std::runtime_error("input row has extra fields");
            fields[index] = std::string_view(line).substr(begin);
        } else {
            if (end == std::string::npos) throw std::runtime_error("input row has too few fields");
            fields[index] = std::string_view(line).substr(begin, end - begin);
            begin = end + 1U;
        }
    }
    return fields;
}

static int parse_positive_int(std::string_view text, const char* label) {
    int value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc() || parsed.ptr != text.data() + text.size() || value <= 0) {
        throw std::runtime_error(std::string("invalid ") + label);
    }
    return value;
}

static size_t parse_positive_size(std::string_view text, const char* label) {
    size_t value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc() || parsed.ptr != text.data() + text.size() || value == 0U) {
        throw std::runtime_error(std::string("invalid ") + label);
    }
    return value;
}

static int parse_nonnegative_int(std::string_view text, const char* label) {
    int value = -1;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc() || parsed.ptr != text.data() + text.size() || value < 0) {
        throw std::runtime_error(std::string("invalid ") + label);
    }
    return value;
}

static double parse_finite_double(std::string_view text, const char* label) {
    // libstdc++ 9 (the server's default compiler) does not provide the
    // floating-point overload of std::from_chars.  Use strtod here while
    // preserving the same full-token and finite-value checks.
    const std::string token(text);
    char* end = nullptr;
    errno = 0;
    const double value = std::strtod(token.c_str(), &end);
    if (end != token.c_str() + token.size() || errno == ERANGE || !std::isfinite(value)) {
        throw std::runtime_error(std::string("invalid ") + label);
    }
    return value;
}

static bool gz_read_line(gzFile handle, std::string& output) {
    output.clear();
    std::array<char, 65536> buffer{};
    while (true) {
        char* observed = gzgets(handle, buffer.data(), static_cast<int>(buffer.size()));
        if (observed == nullptr) {
            if (output.empty()) return false;
            int error_number = Z_OK;
            const char* message = gzerror(handle, &error_number);
            if (error_number != Z_OK && error_number != Z_STREAM_END) {
                throw std::runtime_error(std::string("gzip read failed: ") + message);
            }
            break;
        }
        output.append(observed);
        if (!output.empty() && output.back() == '\n') break;
    }
    if (!output.empty() && output.back() == '\n') output.pop_back();
    if (!output.empty() && output.back() == '\r') output.pop_back();
    return true;
}

static OrderRecord parse_order_record(const std::string& line) {
    static const std::string prefix = "{\"cluster_id\":\"";
    static const std::string entropy_marker = "\",\"entropy_bits\":\"";
    static const std::string permutation_marker = "\",\"permutation\":[";
    if (line.rfind(prefix, 0) != 0 || line.size() < prefix.size() + 4U || line.back() != '}') {
        throw std::runtime_error("invalid order JSON envelope");
    }
    const size_t cluster_end = line.find(entropy_marker, prefix.size());
    if (cluster_end == std::string::npos) throw std::runtime_error("order JSON missing entropy marker");
    const size_t permutation_begin_marker = line.find(permutation_marker, cluster_end + entropy_marker.size());
    if (permutation_begin_marker == std::string::npos) {
        throw std::runtime_error("order JSON missing permutation marker");
    }
    const size_t array_begin = permutation_begin_marker + permutation_marker.size();
    if (line.size() < array_begin + 2U || line[line.size() - 2U] != ']') {
        throw std::runtime_error("order JSON invalid permutation terminator");
    }
    OrderRecord record;
    record.cluster_id = line.substr(prefix.size(), cluster_end - prefix.size());
    const char* cursor = line.data() + array_begin;
    const char* end = line.data() + line.size() - 2U;
    if (cursor == end) return record;
    while (cursor < end) {
        size_t value = 0;
        const auto parsed = std::from_chars(cursor, end, value);
        if (parsed.ec != std::errc() || parsed.ptr == cursor) {
            throw std::runtime_error("invalid permutation integer");
        }
        record.permutation.push_back(value);
        cursor = parsed.ptr;
        if (cursor == end) break;
        if (*cursor != ',') throw std::runtime_error("invalid permutation delimiter");
        ++cursor;
    }
    return record;
}

static std::vector<ClusterTask> load_tasks(const fs::path& reads_path, const fs::path& order_path,
                                           const std::string& dataset, int minimum_design_length,
                                           size_t& total_reads) {
    std::ifstream reads(reads_path);
    if (!reads) throw std::runtime_error("cannot open reads input");
    static std::vector<char> reads_buffer(8U * 1024U * 1024U);
    reads.rdbuf()->pubsetbuf(reads_buffer.data(), reads_buffer.size());
    std::string line;
    if (!std::getline(reads, line) ||
        line != "dataset_id\tcluster_id\tread_id\tread_sequence\toriginal_read_sequence\tread_quality\tdesign_length") {
        throw std::runtime_error("unexpected reads header");
    }
    GzipReader order(order_path);

    std::vector<ClusterTask> tasks;
    tasks.reserve(110000U);
    std::string current_cluster;
    std::vector<std::string> current_reads;
    std::unordered_set<std::string> current_read_ids;
    std::unordered_set<std::string> observed_clusters;
    int current_design_length = 0;
    std::string order_line;

    auto finalize_cluster = [&]() {
        if (current_cluster.empty()) return;
        if (!gz_read_line(order.get(), order_line)) {
            throw std::runtime_error("order pack ended before reads input");
        }
        OrderRecord record = parse_order_record(order_line);
        if (record.cluster_id != current_cluster) {
            throw std::runtime_error("order/input cluster mismatch: " + record.cluster_id + " != " + current_cluster);
        }
        if (record.permutation.size() != current_reads.size()) {
            throw std::runtime_error("permutation/read count mismatch at cluster " + current_cluster);
        }
        std::vector<unsigned char> seen(current_reads.size(), 0U);
        std::vector<std::string> ordered;
        ordered.reserve(current_reads.size());
        for (size_t index : record.permutation) {
            if (index >= current_reads.size() || seen[index] != 0U) {
                throw std::runtime_error("invalid permutation at cluster " + current_cluster);
            }
            seen[index] = 1U;
            ordered.push_back(std::move(current_reads[index]));
        }
        total_reads += ordered.size();
        tasks.push_back({current_cluster, std::move(ordered), current_design_length});
        current_cluster.clear();
        current_reads.clear();
        current_read_ids.clear();
        current_design_length = 0;
    };

    size_t input_rows = 0;
    while (std::getline(reads, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto fields = split_input_row(line);
        if (fields[0] != dataset) throw std::runtime_error("unexpected input dataset");
        if (fields[0].empty() || fields[1].empty() || fields[2].empty() || fields[3].empty()) {
            throw std::runtime_error("empty required input field");
        }
        const std::string cluster(fields[1]);
        const int design_length = parse_positive_int(fields[6], "design_length");
        if (design_length < minimum_design_length || design_length > 1000000) {
            throw std::runtime_error("design_length outside safe range");
        }
        if (current_cluster.empty()) {
            if (!observed_clusters.insert(cluster).second) {
                throw std::runtime_error("duplicate cluster id " + cluster);
            }
            current_cluster = cluster;
            current_design_length = design_length;
        } else if (cluster != current_cluster) {
            finalize_cluster();
            if (!observed_clusters.insert(cluster).second) {
                throw std::runtime_error("duplicate cluster id " + cluster);
            }
            current_cluster = cluster;
            current_design_length = design_length;
        } else if (design_length != current_design_length) {
            throw std::runtime_error("design_length changed within cluster " + cluster);
        }
        const std::string read_id(fields[2]);
        if (!current_read_ids.insert(read_id).second) {
            throw std::runtime_error("duplicate read id in cluster " + cluster);
        }
        current_reads.emplace_back(fields[3]);
        ++input_rows;
    }
    if (!reads.eof()) throw std::runtime_error("reads input ended without clean EOF");
    finalize_cluster();
    if (input_rows == 0 || tasks.empty()) throw std::runtime_error("empty reads input");
    if (gz_read_line(order.get(), order_line)) throw std::runtime_error("order pack has extra rows");
    int gzip_error = Z_OK;
    const char* gzip_message = gzerror(order.get(), &gzip_error);
    if (gzip_error != Z_OK && gzip_error != Z_STREAM_END) {
        throw std::runtime_error(std::string("gzip close audit failed: ") + gzip_message);
    }
    order.close_checked();
    return tasks;
}

static int config_k(const std::string& config_name) {
    constexpr std::string_view prefix = "kmc_k";
    if (config_name.rfind(prefix, 0) != 0) {
        throw std::runtime_error("invalid config prefix");
    }
    const size_t end = config_name.find("_b", prefix.size());
    if (end == std::string::npos) throw std::runtime_error("invalid config k");
    int value = 0;
    const std::string_view text(config_name.data() + prefix.size(), end - prefix.size());
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc() || parsed.ptr != text.data() + text.size() ||
        value < 5 || value > 15) {
        throw std::runtime_error("unsupported config k");
    }
    return value;
}

static ClusterResult process_cluster(
    const ClusterTask& task, WorkerContext& context, size_t worker_id,
    bool profile_enabled, bool hot_cache_enabled = false) {
    const auto cluster_started = std::chrono::steady_clock::now();
    auto profile_started = cluster_started;
    double model_reset_seconds = 0.0;
    if (profile_enabled) profile_started = std::chrono::steady_clock::now();
    context.model.reset(task.design_length);
    // The ordinary cluster path reuses one kernel per worker.  Reserve the
    // cluster's metadata capacity before the first read so read multiplicity
    // hashing and length samples do not repeatedly rehash/grow while the
    // evidence ledger is being accumulated.
    context.model.reserve_evidence(task.reads.size());
    if (hot_cache_enabled && !context.hot_cache_enabled) {
        context.model.enable_hot_cache();
        context.hot_cache_enabled = true;
    }
    if (profile_enabled) {
        model_reset_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - profile_started).count();
    }
    double canonical_seconds = 0.0;
    double add_read_seconds = 0.0;

    ClusterResult output;
    output.cluster_id = task.cluster_id;
    output.total_reads = task.reads.size();
    output.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
    output.worker_id = worker_id;
    output.model_reset_seconds = model_reset_seconds;

    std::string sequence;
    for (size_t index = 0; index < task.reads.size(); ++index) {
        const bool is_last = index + 1U == task.reads.size();
        if (profile_enabled) profile_started = std::chrono::steady_clock::now();
        acor_canonical_into(task.reads[index], sequence);
        if (profile_enabled) {
            const auto finished = std::chrono::steady_clock::now();
            canonical_seconds += std::chrono::duration<double>(finished - profile_started).count();
            profile_started = finished;
        }
        const AcorIncrementalResult incremental = context.model.add_read(sequence, is_last);
        if (profile_enabled) {
            const auto finished = std::chrono::steady_clock::now();
            add_read_seconds += std::chrono::duration<double>(finished - profile_started).count();
            profile_started = finished;
        }
        if (!is_last) continue;
        output.candidate = incremental.ranked.empty()
            ? std::string() : incremental.ranked.front().sequence;
        output.stop_index = index + 1U;
        break;
    }
    if (output.stop_index == 0U) {
        throw std::runtime_error("cluster did not finalize: " + task.cluster_id);
    }
    output.completed_empty = output.candidate.empty();
    output.canonical_seconds = canonical_seconds;
    output.add_read_seconds = add_read_seconds;
    output.graph_update_seconds = context.model.graph_update_seconds();
    output.beam_decode_seconds = context.model.beam_decode_seconds();
    output.cluster_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - cluster_started).count();
    return output;
}

static void process_ledger_block(
    const ClusterTask& task, const LedgerBlockTask& block, WorkerContext& context,
    LedgerBlockResult& output, size_t worker_id, bool profile_enabled,
    bool hot_cache_enabled) {
    const auto block_started = std::chrono::steady_clock::now();
    context.model.reset(task.design_length);
    if (hot_cache_enabled) context.model.enable_hot_cache();
    context.model.reserve_evidence(block.end - block.begin);
    double canonical_seconds = 0.0;
    double evidence_seconds = 0.0;
    std::string sequence;
    for (size_t index = block.begin; index < block.end; ++index) {
        auto stage_started = std::chrono::steady_clock::now();
        acor_canonical_into(task.reads[index], sequence);
        if (profile_enabled) {
            const auto finished = std::chrono::steady_clock::now();
            canonical_seconds += std::chrono::duration<double>(finished - stage_started).count();
            stage_started = finished;
        }
        context.model.add_read_evidence_only(sequence);
        if (profile_enabled) {
            evidence_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - stage_started).count();
        }
    }
    output.cluster_index = block.cluster_index;
    output.block_index = block.block_index;
    output.worker_id = worker_id;
    output.reads = block.end - block.begin;
    output.canonical_seconds = canonical_seconds;
    output.evidence_seconds = evidence_seconds;
    output.block_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - block_started).count();
}

static ClusterResult process_ledger_cluster(
    const ClusterTask& task, size_t workers, const std::vector<int>& worker_cpus,
    bool profile_enabled, std::vector<WorkerStats>& worker_stats,
    bool range_merge, bool hot_cache_enabled,
    bool reuse_phase_threads = false, bool sparse_range_merge = false) {
    const auto cluster_started = std::chrono::steady_clock::now();
    const auto ledger_accumulation_started = cluster_started;
    std::vector<std::unique_ptr<WorkerContext>> contexts;
    contexts.reserve(workers);
    for (size_t worker = 0; worker < workers; ++worker) {
        contexts.push_back(std::make_unique<WorkerContext>());
    }
    // Parallel range reduction writes one shared dense array by disjoint
    // address ranges.  The high-k sparse representation instead uses
    // independently allocated hash tables, so reduce those ledgers through
    // the existing deterministic NUMA-aware tree.  Evidence accumulation
    // remains parallel; only the merge implementation changes.
    if (!contexts.empty() && contexts.front()->model.uses_sparse_ledger()) {
        range_merge = false;
        sparse_range_merge = false;
        reuse_phase_threads = false;
    }
    std::vector<LedgerBlockResult> block_results(workers);
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    std::atomic<bool> worker_failed{false};
    auto run_worker = [&](size_t worker) {
        try {
            pin_current_thread(worker_cpus[worker]);
            const size_t begin = task.reads.size() * worker / workers;
            const size_t end = task.reads.size() * (worker + 1U) / workers;
            const LedgerBlockTask block{0U, worker, begin, end};
            process_ledger_block(
                task, block, *contexts[worker], block_results[worker], worker,
                profile_enabled, hot_cache_enabled);
            worker_stats[worker].clusters += end > begin ? 1U : 0U;
            worker_stats[worker].reads_consumed += end - begin;
            worker_stats[worker].cluster_seconds_sum += block_results[worker].block_seconds;
        } catch (...) {
            worker_failed.store(true, std::memory_order_release);
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
        }
    };

    // Experimental phase-reuse path.  The evidence algorithm is unchanged:
    // each worker still owns a disjoint read range and the root ledger still
    // receives deterministic disjoint counter ranges.  The only change is
    // that the same worker threads remain alive across accumulation and range
    // merge, eliminating a second create/affinity/join cycle per large
    // cluster.
    if (reuse_phase_threads && range_merge && workers > 1U) {
        std::vector<AcorKernel> ledgers;
        ledgers.reserve(workers);
        std::mutex phase_mutex;
        std::condition_variable phase_cv;
        size_t accumulation_done = 0U;
        size_t merge_done = 0U;
        bool merge_phase = false;
        bool stop_phase = false;

        auto get_worker_error = [&]() -> std::exception_ptr {
            std::lock_guard<std::mutex> lock(error_mutex);
            return worker_error;
        };
        auto record_phase_error = [&](std::exception_ptr error) {
            worker_failed.store(true, std::memory_order_release);
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::move(error);
        };
        std::vector<std::thread> phase_threads;
        phase_threads.reserve(workers - 1U);
        auto phase_worker = [&](size_t worker) {
            run_worker(worker);
            {
                std::unique_lock<std::mutex> lock(phase_mutex);
                ++accumulation_done;
                phase_cv.notify_all();
                phase_cv.wait(lock, [&]() { return merge_phase || stop_phase; });
                if (stop_phase) return;
            }
            try {
                pin_current_thread(worker_cpus[worker]);
                for (size_t source = 1U; source < workers; ++source) {
                    if (sparse_range_merge) {
                        ledgers[0].merge_counter_sparse_range_from(
                            ledgers[source], worker, workers);
                    } else {
                        ledgers[0].merge_counter_range_from(
                            ledgers[source], worker, workers);
                    }
                }
            } catch (...) {
                record_phase_error(std::current_exception());
            }
            {
                std::lock_guard<std::mutex> lock(phase_mutex);
                ++merge_done;
                phase_cv.notify_all();
            }
        };
        for (size_t worker = 1U; worker < workers; ++worker) {
            phase_threads.emplace_back([&, worker]() { phase_worker(worker); });
        }
        run_worker(0U);
        {
            std::unique_lock<std::mutex> lock(phase_mutex);
            phase_cv.wait(lock, [&]() {
                return accumulation_done == workers - 1U;
            });
        }
        if (worker_failed.load(std::memory_order_acquire)) {
            {
                std::lock_guard<std::mutex> lock(phase_mutex);
                stop_phase = true;
                phase_cv.notify_all();
            }
            for (std::thread& thread : phase_threads) thread.join();
            if (auto error = get_worker_error()) std::rethrow_exception(error);
            throw std::runtime_error("persistent ledger accumulation failed");
        }
        if (hot_cache_enabled) {
            for (auto& context : contexts) context->model.flush_pending_hot_cache();
        }
        const auto ledger_accumulation_finished = std::chrono::steady_clock::now();

        ClusterResult result;
        result.cluster_id = task.cluster_id;
        result.total_reads = task.reads.size();
        result.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
        result.ledger_accumulation_seconds = std::chrono::duration<double>(
            ledger_accumulation_finished - ledger_accumulation_started).count();
        for (size_t worker = 0U; worker < workers; ++worker) {
            ledgers.emplace_back(std::move(contexts[worker]->model));
            result.canonical_seconds += block_results[worker].canonical_seconds;
            result.evidence_ed_seconds += block_results[worker].evidence_seconds;
            result.graph_update_seconds += ledgers.back().graph_update_seconds();
        }

        const auto ledger_merge_started = std::chrono::steady_clock::now();
        {
            std::lock_guard<std::mutex> lock(phase_mutex);
            merge_phase = true;
            phase_cv.notify_all();
        }
        try {
            pin_current_thread(worker_cpus[0]);
            for (size_t source = 1U; source < workers; ++source) {
                if (sparse_range_merge) {
                    ledgers[0].merge_counter_sparse_range_from(
                        ledgers[source], 0U, workers);
                } else {
                    ledgers[0].merge_counter_range_from(ledgers[source], 0U, workers);
                }
            }
        } catch (...) {
            record_phase_error(std::current_exception());
        }
        {
            std::unique_lock<std::mutex> lock(phase_mutex);
            phase_cv.wait(lock, [&]() {
                return merge_done == workers - 1U ||
                    worker_failed.load(std::memory_order_acquire);
            });
            if (worker_failed.load(std::memory_order_acquire)) stop_phase = true;
            phase_cv.notify_all();
        }
        for (std::thread& thread : phase_threads) thread.join();
        if (worker_failed.load(std::memory_order_acquire)) {
            if (auto error = get_worker_error()) std::rethrow_exception(error);
            throw std::runtime_error("persistent ledger merge failed");
        }
        ledgers[0].finalize_counter_range_merge();
        for (size_t source = 1U; source < workers; ++source) {
            ledgers[0].merge_evidence_metadata_from(ledgers[source]);
        }
        const auto ledger_merge_finished = std::chrono::steady_clock::now();
        result.ledger_merge_seconds = std::chrono::duration<double>(
            ledger_merge_finished - ledger_merge_started).count();
        AcorKernel merged = std::move(ledgers[0]);
        const auto ledger_final_decode_started = std::chrono::steady_clock::now();
        const AcorIncrementalResult decoded = merged.decode_final();
        const auto ledger_final_decode_finished = std::chrono::steady_clock::now();
        result.ledger_final_decode_seconds = std::chrono::duration<double>(
            ledger_final_decode_finished - ledger_final_decode_started).count();
        result.candidate = decoded.ranked.empty() ? std::string() : decoded.ranked.front().sequence;
        result.stop_index = task.reads.size();
        result.completed_empty = result.candidate.empty();
        result.beam_decode_seconds = decoded.final_decode_seconds;
        result.cluster_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - cluster_started).count();
        return result;
    }

    std::vector<std::thread> threads;
    threads.reserve(workers > 0U ? workers - 1U : 0U);
    ThreadJoiner thread_joiner(threads);
    for (size_t worker = 1; worker < workers; ++worker) {
        threads.emplace_back([&, worker]() { run_worker(worker); });
    }
    run_worker(0U);
    for (std::thread& thread : threads) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);
    if (hot_cache_enabled) {
        for (auto& context : contexts) context->model.flush_pending_hot_cache();
    }
    const auto ledger_accumulation_finished = std::chrono::steady_clock::now();

    ClusterResult result;
    result.cluster_id = task.cluster_id;
    result.total_reads = task.reads.size();
    result.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
    result.ledger_accumulation_seconds = std::chrono::duration<double>(
        ledger_accumulation_finished - ledger_accumulation_started).count();
    const auto ledger_merge_started = ledger_accumulation_finished;
    std::vector<AcorKernel> ledgers;
    ledgers.reserve(workers);
    for (size_t worker = 0; worker < workers; ++worker) {
        ledgers.emplace_back(std::move(contexts[worker]->model));
        result.canonical_seconds += block_results[worker].canonical_seconds;
        result.evidence_ed_seconds += block_results[worker].evidence_seconds;
        result.graph_update_seconds += ledgers.back().graph_update_seconds();
    }
    if (range_merge) {
        // Each merge thread owns a disjoint range in every counter array.  It
        // reads all worker ledgers but writes only that range of the root,
        // eliminating intermediate full-ledger writes from the tree merge.
        std::vector<std::thread> range_threads;
        range_threads.reserve(workers > 0U ? workers - 1U : 0U);
        for (size_t shard = 1; shard < workers; ++shard) {
            range_threads.emplace_back([&, shard]() {
                pin_current_thread(worker_cpus[shard]);
                for (size_t source = 1; source < workers; ++source) {
                    if (sparse_range_merge) {
                        ledgers[0].merge_counter_sparse_range_from(
                            ledgers[source], shard, workers);
                    } else {
                        ledgers[0].merge_counter_range_from(
                            ledgers[source], shard, workers);
                    }
                }
            });
        }
        for (size_t source = 1; source < workers; ++source) {
            pin_current_thread(worker_cpus[0]);
            if (sparse_range_merge) {
                ledgers[0].merge_counter_sparse_range_from(
                    ledgers[source], 0U, workers);
            } else {
                ledgers[0].merge_counter_range_from(ledgers[source], 0U, workers);
            }
        }
        for (std::thread& thread : range_threads) thread.join();
        ledgers[0].finalize_counter_range_merge();
        for (size_t source = 1; source < workers; ++source) {
            ledgers[0].merge_evidence_metadata_from(ledgers[source]);
        }
    } else {
    // Merge local ledgers with NUMA locality: reduce each socket first, then
    // perform the small number of cross-socket merges.  Pairing remains fixed
    // and deterministic; the left ledger owns each merge result.
    auto merge_group = [&](std::vector<size_t> active) -> size_t {
        size_t merge_round = 0;
        while (active.size() > 1U) {
        std::vector<std::thread> merge_threads;
        merge_threads.reserve(active.size() / 2U);
        const size_t pair_count = active.size() / 2U;
        for (size_t pair = 0; pair < pair_count; ++pair) {
            const size_t left = active[2U * pair];
            const size_t right = active[2U * pair + 1U];
            merge_threads.emplace_back([&, left, right, pair, merge_round]() {
                pin_current_thread(worker_cpus[left]);
                ledgers[left].merge_evidence_from(ledgers[right]);
            });
        }
        for (std::thread& thread : merge_threads) thread.join();
        std::vector<size_t> next_active;
        next_active.reserve((active.size() + 1U) / 2U);
        for (size_t pair = 0; pair < pair_count; ++pair) next_active.push_back(active[2U * pair]);
        if (active.size() % 2U != 0U) next_active.push_back(active.back());
        active.swap(next_active);
        ++merge_round;
        }
        return active.front();
    };
    std::vector<std::vector<size_t>> socket_groups;
    for (size_t worker = 0; worker < workers; ++worker) {
        const int socket = cpu_location(worker_cpus[worker]).socket;
        auto group = std::find_if(socket_groups.begin(), socket_groups.end(),
            [socket, &worker_cpus](const std::vector<size_t>& members) {
                return !members.empty() &&
                    cpu_location(worker_cpus[members.front()]).socket == socket;
            });
        if (group == socket_groups.end()) socket_groups.push_back({worker});
        else group->push_back(worker);
    }
    std::vector<size_t> socket_roots;
    socket_roots.reserve(socket_groups.size());
    for (const auto& group : socket_groups) socket_roots.push_back(merge_group(group));
    const size_t merged_root = socket_roots.size() > 1U
        ? merge_group(socket_roots) : socket_roots.front();
    if (merged_root != 0U) {
        ledgers[0] = std::move(ledgers[merged_root]);
    }
    }
    const auto ledger_merge_finished = std::chrono::steady_clock::now();
    result.ledger_merge_seconds = std::chrono::duration<double>(
        ledger_merge_finished - ledger_merge_started).count();
    AcorKernel merged = std::move(ledgers[0]);
    const auto ledger_final_decode_started = std::chrono::steady_clock::now();
    const AcorIncrementalResult decoded = merged.decode_final();
    const auto ledger_final_decode_finished = std::chrono::steady_clock::now();
    result.ledger_final_decode_seconds = std::chrono::duration<double>(
        ledger_final_decode_finished - ledger_final_decode_started).count();
    result.candidate = decoded.ranked.empty() ? std::string() : decoded.ranked.front().sequence;
    result.stop_index = task.reads.size();
    result.completed_empty = result.candidate.empty();
    result.beam_decode_seconds = decoded.final_decode_seconds;
    result.cluster_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - cluster_started).count();
    return result;
}

// Persistent block-pool variant for the large-cluster path.  Unlike
// process_ledger_cluster(), this function creates the worker threads once for
// all eligible clusters.  Each worker reuses one WorkerContext while local
// ledgers are moved into block_results.  The final reduction is deterministic:
// blocks are merged in increasing block order and metadata is merged in the
// same order as the input read ranges.
struct PersistentLedgerBlockState {
    LedgerBlockTask task;
    std::unique_ptr<AcorKernel> ledger;
    LedgerBlockResult metrics;
};

static void run_persistent_ledger_block_pool(
    const std::vector<size_t>& large_indices,
    const std::vector<ClusterTask>& tasks,
    std::vector<ClusterResult>& results,
    size_t workers, const std::vector<int>& worker_cpus,
    bool profile_enabled, size_t block_threshold,
    std::vector<WorkerStats>& worker_stats) {
    if (large_indices.empty()) return;

    std::vector<PersistentLedgerBlockState> block_states;
    std::vector<std::vector<size_t>> cluster_blocks(tasks.size());
    for (const size_t cluster_index : large_indices) {
        const size_t read_count = tasks[cluster_index].reads.size();
        const size_t target_blocks = std::max<size_t>(
            2U, std::min(workers, (read_count + block_threshold - 1U) / block_threshold));
        const size_t block_count = std::min(target_blocks, read_count);
        auto& indices = cluster_blocks[cluster_index];
        indices.reserve(block_count);
        for (size_t block_index = 0; block_index < block_count; ++block_index) {
            const size_t begin = read_count * block_index / block_count;
            const size_t end = read_count * (block_index + 1U) / block_count;
            if (begin == end) continue;
            indices.push_back(block_states.size());
            PersistentLedgerBlockState state;
            state.task = LedgerBlockTask{cluster_index, block_index, begin, end};
            state.metrics.cluster_index = cluster_index;
            state.metrics.block_index = block_index;
            state.metrics.reads = end - begin;
            block_states.push_back(std::move(state));
        }
    }

    std::atomic<size_t> next{0U};
    std::atomic<bool> failed{false};
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto worker_loop = [&](size_t worker) {
        WorkerContext context;
        try {
            pin_current_thread(worker_cpus[worker]);
            while (!failed.load(std::memory_order_relaxed)) {
                const size_t block_position = next.fetch_add(1U, std::memory_order_relaxed);
                if (block_position >= block_states.size()) break;
                auto& state = block_states[block_position];
                const auto& task = tasks[state.task.cluster_index];
                const auto block_started = std::chrono::steady_clock::now();
                context.model.reset(task.design_length);
                context.model.reserve_evidence(state.task.end - state.task.begin);
                double canonical_seconds = 0.0;
                double evidence_seconds = 0.0;
                std::string sequence;
                for (size_t index = state.task.begin; index < state.task.end; ++index) {
                    auto stage_started = std::chrono::steady_clock::now();
                    acor_canonical_into(task.reads[index], sequence);
                    if (profile_enabled) {
                        const auto finished = std::chrono::steady_clock::now();
                        canonical_seconds += std::chrono::duration<double>(
                            finished - stage_started).count();
                        stage_started = finished;
                    }
                    context.model.add_read_evidence_only(sequence);
                    if (profile_enabled) {
                        evidence_seconds += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - stage_started).count();
                    }
                }
                state.metrics.worker_id = worker;
                state.metrics.canonical_seconds = canonical_seconds;
                state.metrics.evidence_seconds = evidence_seconds;
                state.metrics.block_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - block_started).count();
                state.ledger = std::make_unique<AcorKernel>(std::move(context.model));
                // Moving the model into the block result leaves the worker's
                // context empty. Recreate it before the worker claims the
                // next block; otherwise the next reset() dereferences a null
                // implementation and crashes on multi-block clusters.
                context.model = AcorKernel();
                ++worker_stats[worker].clusters;
                worker_stats[worker].reads_consumed += state.metrics.reads;
                worker_stats[worker].cluster_seconds_sum += state.metrics.block_seconds;
            }
        } catch (...) {
            failed.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(workers > 0U ? workers - 1U : 0U);
    ThreadJoiner thread_joiner(threads);
    for (size_t worker = 1; worker < workers; ++worker) {
        threads.emplace_back([&, worker]() { worker_loop(worker); });
    }
    worker_loop(0U);
    for (std::thread& thread : threads) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);

    // Reduce each eligible cluster independently.  This phase is deliberately
    // deterministic and kept separate from the block pool so that any speedup
    // can be attributed to evidence accumulation rather than altered decoding.
    for (const size_t cluster_index : large_indices) {
        const auto& block_indices = cluster_blocks[cluster_index];
        if (block_indices.empty()) throw std::runtime_error("empty persistent ledger cluster");
        const auto& task = tasks[cluster_index];
        const auto merge_started = std::chrono::steady_clock::now();
        AcorKernel merged = std::move(*block_states[block_indices.front()].ledger);
        double canonical_seconds = 0.0;
        double evidence_seconds = 0.0;
        for (const size_t state_index : block_indices) {
            const auto& state = block_states[state_index];
            canonical_seconds += state.metrics.canonical_seconds;
            evidence_seconds += state.metrics.evidence_seconds;
        }
        for (size_t position = 1; position < block_indices.size(); ++position) {
            const auto& state = block_states[block_indices[position]];
            merged.merge_evidence_from(*state.ledger);
        }
        for (size_t position = 1; position < block_indices.size(); ++position) {
            const auto& state = block_states[block_indices[position]];
            merged.merge_evidence_metadata_from(*state.ledger);
        }
        const auto merge_finished = std::chrono::steady_clock::now();
        const auto decode_started = merge_finished;
        const AcorIncrementalResult decoded = merged.decode_final();
        const auto decode_finished = std::chrono::steady_clock::now();

        ClusterResult result;
        result.cluster_id = task.cluster_id;
        result.total_reads = task.reads.size();
        result.stop_index = task.reads.size();
        result.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
        result.candidate = decoded.ranked.empty() ? std::string() : decoded.ranked.front().sequence;
        result.completed_empty = result.candidate.empty();
        result.canonical_seconds = canonical_seconds;
        result.evidence_ed_seconds = evidence_seconds;
        result.graph_update_seconds = merged.graph_update_seconds();
        result.beam_decode_seconds = decoded.final_decode_seconds;
        result.ledger_accumulation_seconds = 0.0;
        for (const size_t state_index : block_indices) {
            result.ledger_accumulation_seconds += block_states[state_index].metrics.block_seconds;
        }
        result.ledger_merge_seconds = std::chrono::duration<double>(
            merge_finished - merge_started).count();
        result.ledger_final_decode_seconds = std::chrono::duration<double>(
            decode_finished - decode_started).count();
        result.cluster_seconds = std::chrono::duration<double>(
            decode_finished - merge_started).count() + result.ledger_accumulation_seconds;
        results[cluster_index] = std::move(result);
    }
}

// Experimental owner-sharded evidence path.  Workers build sparse additive
// deltas while scanning disjoint read ranges; counter keys are then assigned
// to a single owner and applied to one final dense ledger.  This removes the
// P-way full-ledger merge while preserving the exact integer evidence state.
static ClusterResult process_sharded_cluster(
    const ClusterTask& task, size_t workers, const std::vector<int>& worker_cpus,
    bool /* profile_enabled */, std::vector<WorkerStats>& worker_stats) {
    const auto cluster_started = std::chrono::steady_clock::now();
    const auto accumulation_started = cluster_started;
    std::vector<std::unique_ptr<SparseEvidenceBlock>> blocks(workers);
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto run_worker = [&](size_t worker) {
        try {
            pin_current_thread(worker_cpus[worker]);
            blocks[worker] = std::make_unique<SparseEvidenceBlock>();
            const size_t begin = task.reads.size() * worker / workers;
            const size_t end = task.reads.size() * (worker + 1U) / workers;
            blocks[worker]->read_counts.reserve(end - begin);
            blocks[worker]->lengths.reserve(end - begin);
            const auto block_started = std::chrono::steady_clock::now();
            std::string sequence;
            for (size_t index = begin; index < end; ++index) {
                acor_canonical_into(task.reads[index], sequence);
                blocks[worker]->add_read(sequence, task.design_length, 9);
            }
            worker_stats[worker].clusters += end > begin ? 1U : 0U;
            worker_stats[worker].reads_consumed += end - begin;
            worker_stats[worker].cluster_seconds_sum += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - block_started).count();
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
        }
    };
    std::vector<std::thread> threads;
    threads.reserve(workers > 0U ? workers - 1U : 0U);
    ThreadJoiner thread_joiner(threads);
    for (size_t worker = 1; worker < workers; ++worker) {
        threads.emplace_back([&, worker]() { run_worker(worker); });
    }
    run_worker(0U);
    for (std::thread& thread : threads) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);
    const auto accumulation_finished = std::chrono::steady_clock::now();

    const auto merge_started = accumulation_finished;
    std::vector<std::vector<AcorCounterDelta>> owner_deltas(workers);
    for (size_t worker = 0; worker < workers; ++worker) {
        for (const auto& item : blocks[worker]->counters) {
            const std::uint8_t kind = static_cast<std::uint8_t>(item.first >> 32U);
            const std::uint32_t index = static_cast<std::uint32_t>(item.first & 0xffffffffULL);
            owner_deltas[sparse_owner(item.first, workers)].push_back({kind, index, item.second});
        }
    }

    AcorKernel merged;
    merged.reset(task.design_length);
    merged.reserve_evidence(task.reads.size());
    std::exception_ptr owner_error;
    std::mutex owner_error_mutex;
    auto apply_owner = [&](size_t owner) {
        try {
            pin_current_thread(worker_cpus[owner]);
            merged.apply_counter_deltas(owner_deltas[owner]);
        } catch (...) {
            std::lock_guard<std::mutex> lock(owner_error_mutex);
            if (!owner_error) owner_error = std::current_exception();
        }
    };
    threads.clear();
    threads.reserve(workers > 0U ? workers - 1U : 0U);
    for (size_t owner = 1; owner < workers; ++owner) {
        threads.emplace_back([&, owner]() { apply_owner(owner); });
    }
    apply_owner(0U);
    for (std::thread& thread : threads) thread.join();
    if (owner_error) std::rethrow_exception(owner_error);
    merged.finalize_counter_deltas();

    // Metadata is additive and kept in worker order so that the floating
    // length summaries and tie-breaking behavior remain deterministic.
    for (size_t worker = 0; worker < workers; ++worker) {
        std::vector<std::pair<std::string, std::uint32_t>> read_counts;
        read_counts.reserve(blocks[worker]->read_counts.size());
        for (const auto& item : blocks[worker]->read_counts) {
            read_counts.push_back(item);
        }
        merged.merge_sparse_metadata(
            read_counts, blocks[worker]->lengths, blocks[worker]->reads_arrived,
            blocks[worker]->exact_m_read_count, blocks[worker]->log_length_sum,
            blocks[worker]->log_length_sq_sum,
            blocks[worker]->total_start_observations);
    }
    const auto merge_finished = std::chrono::steady_clock::now();

    ClusterResult result;
    result.cluster_id = task.cluster_id;
    result.total_reads = task.reads.size();
    result.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
    result.ledger_accumulation_seconds = std::chrono::duration<double>(
        accumulation_finished - accumulation_started).count();
    result.ledger_merge_seconds = std::chrono::duration<double>(
        merge_finished - merge_started).count();
    const auto decoded = merged.decode_final();
    result.ledger_final_decode_seconds = decoded.final_decode_seconds;
    result.candidate = decoded.ranked.empty() ? std::string() : decoded.ranked.front().sequence;
    result.stop_index = task.reads.size();
    result.completed_empty = result.candidate.empty();
    result.beam_decode_seconds = decoded.final_decode_seconds;
    result.cluster_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - cluster_started).count();
    return result;
}

// Shared atomic dense-ledger reduction.  This deliberately removes the P-way
// dense-ledger replication and the subsequent full-ledger merge.  It is an
// exact diagnostic route: graph counters use relaxed atomic increments, while
// all non-graph evidence remains private and is merged in worker order.
static ClusterResult process_atomic_shared_cluster(
    const ClusterTask& task, size_t workers, const std::vector<int>& worker_cpus,
    bool profile_enabled, std::vector<WorkerStats>& worker_stats) {
    const auto cluster_started = std::chrono::steady_clock::now();
    const auto accumulation_started = cluster_started;
    std::vector<AtomicMetadataBlock> metadata(workers);
    for (size_t worker = 0; worker < workers; ++worker) {
        const size_t begin = task.reads.size() * worker / workers;
        const size_t end = task.reads.size() * (worker + 1U) / workers;
        metadata[worker].read_counts.reserve(end - begin + 1U);
        metadata[worker].lengths.reserve(end - begin);
    }

    AcorKernel merged;
    merged.reset(task.design_length);
    merged.reserve_evidence(task.reads.size());
    merged.enable_atomic_updates();

    std::vector<double> canonical_seconds(workers, 0.0);
    std::vector<double> evidence_seconds(workers, 0.0);
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto run_worker = [&](size_t worker) {
        try {
            pin_current_thread(worker_cpus[worker]);
            const size_t begin = task.reads.size() * worker / workers;
            const size_t end = task.reads.size() * (worker + 1U) / workers;
            const auto worker_started = std::chrono::steady_clock::now();
            std::string sequence;
            for (size_t index = begin; index < end; ++index) {
                const auto canonical_started = std::chrono::steady_clock::now();
                acor_canonical_into(task.reads[index], sequence);
                const auto canonical_finished = std::chrono::steady_clock::now();
                merged.add_read_graph_atomic(sequence);
                metadata[worker].add_read(sequence, task.design_length, 9);
                const auto finished = std::chrono::steady_clock::now();
                if (profile_enabled) {
                    canonical_seconds[worker] += std::chrono::duration<double>(
                        canonical_finished - canonical_started).count();
                    evidence_seconds[worker] += std::chrono::duration<double>(
                        finished - canonical_finished).count();
                }
            }
            worker_stats[worker].clusters += end > begin ? 1U : 0U;
            worker_stats[worker].reads_consumed += end - begin;
            worker_stats[worker].cluster_seconds_sum += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - worker_started).count();
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
        }
    };
    std::vector<std::thread> threads;
    threads.reserve(workers > 0U ? workers - 1U : 0U);
    ThreadJoiner thread_joiner(threads);
    for (size_t worker = 1; worker < workers; ++worker) {
        threads.emplace_back([&, worker]() { run_worker(worker); });
    }
    run_worker(0U);
    for (std::thread& thread : threads) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);

    const auto accumulation_finished = std::chrono::steady_clock::now();
    const auto merge_started = accumulation_finished;
    merged.flush_atomic_updates();
    for (size_t worker = 0; worker < workers; ++worker) {
        std::vector<std::pair<std::string, std::uint32_t>> read_counts;
        read_counts.reserve(metadata[worker].read_counts.size());
        for (const auto& item : metadata[worker].read_counts) read_counts.push_back(item);
        merged.merge_sparse_metadata(
            read_counts, metadata[worker].lengths, metadata[worker].reads_arrived,
            metadata[worker].exact_m_read_count, metadata[worker].log_length_sum,
            metadata[worker].log_length_sq_sum,
            metadata[worker].total_start_observations);
    }
    const auto merge_finished = std::chrono::steady_clock::now();

    const AcorIncrementalResult decoded = merged.decode_final();
    ClusterResult result;
    result.cluster_id = task.cluster_id;
    result.total_reads = task.reads.size();
    result.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
    result.canonical_seconds = std::accumulate(canonical_seconds.begin(), canonical_seconds.end(), 0.0);
    result.evidence_ed_seconds = std::accumulate(evidence_seconds.begin(), evidence_seconds.end(), 0.0);
    result.ledger_accumulation_seconds = std::chrono::duration<double>(
        accumulation_finished - accumulation_started).count();
    result.ledger_merge_seconds = std::chrono::duration<double>(
        merge_finished - merge_started).count();
    result.ledger_final_decode_seconds = decoded.final_decode_seconds;
    result.candidate = decoded.ranked.empty() ? std::string() : decoded.ranked.front().sequence;
    result.stop_index = task.reads.size();
    result.completed_empty = result.candidate.empty();
    result.beam_decode_seconds = decoded.final_decode_seconds;
    result.cluster_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - cluster_started).count();
    return result;
}

static void tile_descriptor(std::size_t id, std::uint8_t& kind, std::uint32_t& begin) {
    if (id < TileLedger::tiles_per_node) {
        kind = 0U;
        begin = static_cast<std::uint32_t>(id * TileLedger::tile_size);
    } else if (id < TileLedger::tiles_per_node + TileLedger::tiles_per_edge) {
        kind = 1U;
        begin = static_cast<std::uint32_t>(
            (id - TileLedger::tiles_per_node) * TileLedger::tile_size);
    } else if (id < TileLedger::tiles_per_node * 2U + TileLedger::tiles_per_edge) {
        kind = 2U;
        begin = static_cast<std::uint32_t>(
            (id - TileLedger::tiles_per_node - TileLedger::tiles_per_edge) *
            TileLedger::tile_size);
    } else if (id < TileLedger::tiles_per_node * 3U + TileLedger::tiles_per_edge) {
        kind = 3U;
        begin = static_cast<std::uint32_t>(
            (id - TileLedger::tiles_per_node * 2U - TileLedger::tiles_per_edge) *
            TileLedger::tile_size);
    } else {
        kind = 4U;
        begin = static_cast<std::uint32_t>(
            (id - TileLedger::tiles_per_node * 3U - TileLedger::tiles_per_edge) *
            TileLedger::tile_size);
    }
}

// Sparse owner-tile reduction.  Workers update only local fixed-size tile
// payloads.  After accumulation, each tile has one owner thread, so the
// final ledger receives disjoint writes without atomics or a full raw-event
// sort.  This is the first owner-reduction route that does not create one
// object per k-mer event.
static ClusterResult process_owner_tile_cluster(
    const ClusterTask& task, size_t workers, const std::vector<int>& worker_cpus,
    bool profile_enabled, std::vector<WorkerStats>& worker_stats) {
    const auto cluster_started = std::chrono::steady_clock::now();
    const auto accumulation_started = cluster_started;
    std::vector<TileLedger> ledgers(workers);
    std::vector<AtomicMetadataBlock> metadata(workers);
    for (size_t worker = 0; worker < workers; ++worker) {
        const size_t begin = task.reads.size() * worker / workers;
        const size_t end = task.reads.size() * (worker + 1U) / workers;
        metadata[worker].read_counts.reserve(end - begin + 1U);
        metadata[worker].lengths.reserve(end - begin);
    }
    std::vector<double> canonical_seconds(workers, 0.0);
    std::vector<double> evidence_seconds(workers, 0.0);
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto run_worker = [&](size_t worker) {
        try {
            pin_current_thread(worker_cpus[worker]);
            const size_t begin = task.reads.size() * worker / workers;
            const size_t end = task.reads.size() * (worker + 1U) / workers;
            const auto worker_started = std::chrono::steady_clock::now();
            std::string sequence;
            for (size_t index = begin; index < end; ++index) {
                const auto canonical_started = std::chrono::steady_clock::now();
                acor_canonical_into(task.reads[index], sequence);
                const auto canonical_finished = std::chrono::steady_clock::now();
                ledgers[worker].add_sequence(sequence, 9);
                metadata[worker].add_read(sequence, task.design_length, 9);
                const auto finished = std::chrono::steady_clock::now();
                if (profile_enabled) {
                    canonical_seconds[worker] += std::chrono::duration<double>(
                        canonical_finished - canonical_started).count();
                    evidence_seconds[worker] += std::chrono::duration<double>(
                        finished - canonical_finished).count();
                }
            }
            worker_stats[worker].clusters += end > begin ? 1U : 0U;
            worker_stats[worker].reads_consumed += end - begin;
            worker_stats[worker].cluster_seconds_sum += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - worker_started).count();
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
        }
    };
    std::vector<std::thread> threads;
    threads.reserve(workers > 0U ? workers - 1U : 0U);
    ThreadJoiner thread_joiner(threads);
    for (size_t worker = 1; worker < workers; ++worker) {
        threads.emplace_back([&, worker]() { run_worker(worker); });
    }
    run_worker(0U);
    for (std::thread& thread : threads) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);
    const auto accumulation_finished = std::chrono::steady_clock::now();

    const auto merge_started = accumulation_finished;
    AcorKernel merged;
    merged.reset(task.design_length);
    merged.reserve_evidence(task.reads.size());
    std::exception_ptr owner_error;
    std::mutex owner_error_mutex;
    auto reduce_owner = [&](size_t owner) {
        try {
            pin_current_thread(worker_cpus[owner]);
            for (std::size_t tile = owner; tile < TileLedger::total_tiles;
                 tile += workers) {
                std::uint8_t kind = 0U;
                std::uint32_t begin = 0U;
                tile_descriptor(tile, kind, begin);
                for (size_t worker = 0; worker < workers; ++worker) {
                    const auto& payload = ledgers[worker].tiles[tile];
                    if (!payload) continue;
                    merged.apply_counter_tile(
                        kind, begin, payload->counts.data(), TileLedger::tile_size);
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(owner_error_mutex);
            if (!owner_error) owner_error = std::current_exception();
        }
    };
    threads.clear();
    threads.reserve(workers > 0U ? workers - 1U : 0U);
    for (size_t owner = 1; owner < workers; ++owner) {
        threads.emplace_back([&, owner]() { reduce_owner(owner); });
    }
    reduce_owner(0U);
    for (std::thread& thread : threads) thread.join();
    if (owner_error) std::rethrow_exception(owner_error);
    for (size_t worker = 0; worker < workers; ++worker) {
        std::vector<std::pair<std::string, std::uint32_t>> read_counts;
        read_counts.reserve(metadata[worker].read_counts.size());
        for (const auto& item : metadata[worker].read_counts) read_counts.push_back(item);
        merged.merge_sparse_metadata(
            read_counts, metadata[worker].lengths, metadata[worker].reads_arrived,
            metadata[worker].exact_m_read_count, metadata[worker].log_length_sum,
            metadata[worker].log_length_sq_sum,
            metadata[worker].total_start_observations);
    }
    merged.finalize_counter_deltas();
    const auto merge_finished = std::chrono::steady_clock::now();
    const AcorIncrementalResult decoded = merged.decode_final();

    ClusterResult result;
    result.cluster_id = task.cluster_id;
    result.total_reads = task.reads.size();
    result.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
    result.canonical_seconds = std::accumulate(canonical_seconds.begin(), canonical_seconds.end(), 0.0);
    result.evidence_ed_seconds = std::accumulate(evidence_seconds.begin(), evidence_seconds.end(), 0.0);
    result.ledger_accumulation_seconds = std::chrono::duration<double>(
        accumulation_finished - accumulation_started).count();
    result.ledger_merge_seconds = std::chrono::duration<double>(
        merge_finished - merge_started).count();
    result.ledger_final_decode_seconds = decoded.final_decode_seconds;
    result.candidate = decoded.ranked.empty() ? std::string() : decoded.ranked.front().sequence;
    result.stop_index = task.reads.size();
    result.completed_empty = result.candidate.empty();
    result.beam_decode_seconds = decoded.final_decode_seconds;
    result.cluster_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - cluster_started).count();
    return result;
}

// Batched owner-sharded reduction.  Workers scan disjoint read ranges into
// bounded append-only update buffers.  After each batch, updates are routed
// to a deterministic owner, sorted by (counter kind, index), run-length
// reduced, and applied to disjoint counter ownership domains in one final
// ledger.  No worker ever owns a complete dense ledger.
static ClusterResult process_batch_owner_cluster(
    const ClusterTask& task, size_t workers, const std::vector<int>& worker_cpus,
    bool profile_enabled, std::vector<WorkerStats>& worker_stats) {
    constexpr size_t batch_reads_per_worker = 4096U;
    const auto cluster_started = std::chrono::steady_clock::now();
    const auto accumulation_started = cluster_started;
    std::vector<BatchEvidenceBlock> blocks(workers);
    for (size_t worker = 0; worker < workers; ++worker) {
        blocks[worker].read_counts.reserve(task.reads.size() / workers + 1U);
        blocks[worker].lengths.reserve(task.reads.size() / workers + 1U);
        blocks[worker].deltas.reserve(batch_reads_per_worker * 900U);
    }
    AcorKernel emitter;
    emitter.reset(task.design_length);
    AcorKernel merged;
    merged.reset(task.design_length);
    merged.reserve_evidence(task.reads.size());

    double canonical_seconds = 0.0;
    double evidence_seconds = 0.0;
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    for (size_t batch_begin = 0; batch_begin < task.reads.size();
         batch_begin += batch_reads_per_worker * workers) {
        const size_t batch_end = std::min(
            task.reads.size(), batch_begin + batch_reads_per_worker * workers);
        for (BatchEvidenceBlock& block : blocks) block.clear_batch();
        auto run_worker = [&](size_t worker) {
            try {
                pin_current_thread(worker_cpus[worker]);
                const size_t begin = std::min(
                    batch_end, batch_begin + batch_reads_per_worker * worker);
            const size_t end = std::min(
                    batch_end, begin + batch_reads_per_worker);
                std::string sequence;
                for (size_t index = begin; index < end; ++index) {
                    const auto canonical_started = std::chrono::steady_clock::now();
                    acor_canonical_into(task.reads[index], sequence);
                    if (profile_enabled) {
                        blocks[worker].canonical_seconds += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - canonical_started).count();
                    }
                    blocks[worker].add_read(
                        sequence, task.design_length, emitter, profile_enabled);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };
        std::vector<std::thread> threads;
        threads.reserve(workers > 0U ? workers - 1U : 0U);
        ThreadJoiner thread_joiner(threads);
        for (size_t worker = 1; worker < workers; ++worker) {
            threads.emplace_back([&, worker]() { run_worker(worker); });
        }
        run_worker(0U);
        for (std::thread& thread : threads) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        const auto accumulation_batch_finished = std::chrono::steady_clock::now();
        std::vector<std::vector<AcorCounterDelta>> owner_updates(workers);
        for (size_t worker = 0; worker < workers; ++worker) {
            for (const AcorCounterDelta& delta : blocks[worker].deltas) {
                const std::uint64_t key =
                    (static_cast<std::uint64_t>(delta.kind) << 32U) | delta.index;
                owner_updates[sparse_owner(key, workers)].push_back(delta);
            }
        }

        std::exception_ptr owner_error;
        std::mutex owner_error_mutex;
        auto reduce_owner = [&](size_t owner) {
            try {
                pin_current_thread(worker_cpus[owner]);
                std::vector<AcorCounterDelta>& updates = owner_updates[owner];
                std::sort(updates.begin(), updates.end(), [](
                    const AcorCounterDelta& lhs, const AcorCounterDelta& rhs) {
                    if (lhs.kind != rhs.kind) return lhs.kind < rhs.kind;
                    return lhs.index < rhs.index;
                });
                std::vector<AcorCounterDelta> reduced;
                reduced.reserve(updates.size());
                for (const AcorCounterDelta& delta : updates) {
                    if (!reduced.empty() && reduced.back().kind == delta.kind &&
                        reduced.back().index == delta.index) {
                        if (reduced.back().amount >
                            std::numeric_limits<std::uint32_t>::max() - delta.amount) {
                            throw std::runtime_error("batched counter reduction overflow");
                        }
                        reduced.back().amount += delta.amount;
                    } else {
                        reduced.push_back(delta);
                    }
                }
                merged.apply_counter_deltas(reduced);
            } catch (...) {
                std::lock_guard<std::mutex> lock(owner_error_mutex);
                if (!owner_error) owner_error = std::current_exception();
            }
        };
        threads.clear();
        threads.reserve(workers > 0U ? workers - 1U : 0U);
        for (size_t owner = 1; owner < workers; ++owner) {
            threads.emplace_back([&, owner]() { reduce_owner(owner); });
        }
        reduce_owner(0U);
        for (std::thread& thread : threads) thread.join();
        if (owner_error) std::rethrow_exception(owner_error);

        canonical_seconds += std::accumulate(
            blocks.begin(), blocks.end(), 0.0,
            [](double total, const BatchEvidenceBlock& block) {
                return total + block.canonical_seconds;
            });
        evidence_seconds += std::accumulate(
            blocks.begin(), blocks.end(), 0.0,
            [](double total, const BatchEvidenceBlock& block) {
                return total + block.evidence_seconds;
            });
        (void)accumulation_batch_finished;
    }
    const auto accumulation_finished = std::chrono::steady_clock::now();

    const auto merge_started = accumulation_finished;
    for (size_t worker = 0; worker < workers; ++worker) {
        std::vector<std::pair<std::string, std::uint32_t>> read_counts;
        read_counts.reserve(blocks[worker].read_counts.size());
        for (const auto& item : blocks[worker].read_counts) read_counts.push_back(item);
        merged.merge_sparse_metadata(
            read_counts, blocks[worker].lengths, blocks[worker].reads_arrived,
            blocks[worker].exact_m_read_count, blocks[worker].log_length_sum,
            blocks[worker].log_length_sq_sum, blocks[worker].total_start_observations);
    }
    merged.finalize_counter_deltas();
    const auto merge_finished = std::chrono::steady_clock::now();

    for (size_t worker = 0; worker < workers; ++worker) {
        worker_stats[worker].clusters +=
            (task.reads.size() * worker / workers < task.reads.size() * (worker + 1U) / workers)
                ? 1U : 0U;
        worker_stats[worker].reads_consumed +=
            task.reads.size() * (worker + 1U) / workers - task.reads.size() * worker / workers;
    }
    ClusterResult result;
    result.cluster_id = task.cluster_id;
    result.total_reads = task.reads.size();
    result.legacy_prefix_alias = std::min<size_t>(20U, task.reads.size());
    result.canonical_seconds = canonical_seconds;
    result.evidence_ed_seconds = evidence_seconds;
    result.ledger_accumulation_seconds = std::chrono::duration<double>(
        accumulation_finished - accumulation_started).count();
    result.ledger_merge_seconds = std::chrono::duration<double>(
        merge_finished - merge_started).count();
    const auto decoded = merged.decode_final();
    result.ledger_final_decode_seconds = decoded.final_decode_seconds;
    result.candidate = decoded.ranked.empty() ? std::string() : decoded.ranked.front().sequence;
    result.stop_index = task.reads.size();
    result.completed_empty = result.candidate.empty();
    result.beam_decode_seconds = decoded.final_decode_seconds;
    result.cluster_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - cluster_started).count();
    return result;
}

static void write_cpu_audit(const fs::path& path, int control_cpu,
                            const std::vector<int>& worker_cpus, bool allow_hyperthreads) {
    std::ofstream output(path);
    output << "role\tworker_id\tlogical_cpu\tsocket\tphysical_core\tphysical_core_status\n";
    const CpuLocation control = cpu_location(control_cpu);
    output << "control\tNA\t" << control.cpu << '\t' << control.socket << '\t'
           << control.core << "\tCONTROL_SEQUENTIAL_ONLY\n";
    std::set<int> logical_seen;
    std::set<std::pair<int, int>> seen;
    for (size_t index = 0; index < worker_cpus.size(); ++index) {
        const CpuLocation location = cpu_location(worker_cpus[index]);
        const bool logical_unique = logical_seen.insert(location.cpu).second;
        const bool unique = seen.insert({location.socket, location.core}).second;
        if (!logical_unique) throw std::runtime_error("duplicate worker logical cpu");
        output << "worker\t" << index << '\t' << location.cpu << '\t' << location.socket
               << '\t' << location.core << '\t'
               << (unique ? "PASS_UNIQUE" : (allow_hyperthreads ? "HT_SHARED_EXPECTED" : "FAIL_DUPLICATE"))
               << '\n';
        if (!unique && !allow_hyperthreads) throw std::runtime_error("worker hyperthread/core conflict");
    }
    output.flush();
    if (!output) throw std::runtime_error("cpu audit write failed");
    output.close();
    if (!output) throw std::runtime_error("cpu audit close failed");
}

int main(int argc, char** argv) {
    try {
        if (argc != 15) {
            std::cerr << "usage: runner READS.tsv ORDER.jsonl.gz OUTPUT_DIR ROUND_ID INPUT_DATASET OUTPUT_DATASET_ID "
                         "WORKERS CPU_LIST CONTROL_CPU CONFIG MODE MU_STAR DELTA_FLIP {PROFILE|NO_PROFILE}\n";
            return 2;
        }
        const fs::path reads_path = argv[1];
        const fs::path order_path = argv[2];
        const fs::path output_dir = argv[3];
        const std::string round_id = argv[4];
        const std::string dataset = argv[5];
        const std::string output_dataset_id = argv[6];
        const size_t workers = parse_positive_size(argv[7], "workers");
        const std::vector<int> worker_cpus = parse_cpu_list(argv[8]);
        const int control_cpu = parse_nonnegative_int(argv[9], "control cpu");
        if (control_cpu >= CPU_SETSIZE) throw std::runtime_error("control cpu outside CPU_SETSIZE");
        const std::string config_name = argv[10];
        const std::string mode = argv[11];
        const double mu_star = parse_finite_double(argv[12], "mu_star");
        const double delta_flip = parse_finite_double(argv[13], "delta_flip");
        const std::string profile_mode = argv[14];
        const bool profile_enabled = profile_mode == "PROFILE";
        if (!profile_enabled && profile_mode != "NO_PROFILE") throw std::runtime_error("invalid profile mode");
        if (workers == 0U || worker_cpus.size() != workers) {
            throw std::runtime_error("worker/cpu-list size mismatch");
        }
        if (std::find(worker_cpus.begin(), worker_cpus.end(), control_cpu) != worker_cpus.end()) {
            throw std::runtime_error("control cpu must differ from worker cpus");
        }
        const char* configured = std::getenv("ACOR_CONFIG");
        if (configured == nullptr) configured = std::getenv("ACOR_RESEARCH_CONFIG");
        const std::string effective_config = configured == nullptr
            ? "kmc_k9_b16_lognormal" : std::string(configured);
        if (config_name != effective_config) {
            throw std::runtime_error("runner/config environment mismatch");
        }
        const bool sharded_mode = mode == "SHARDED" || mode.rfind("SHARDED_", 0) == 0;
        const bool hot_cache_mode = mode == "LEDGER_HOTCACHE";
        const bool persistent_shardmerge_mode = mode == "LEDGER_PERSIST_SHARDMERGE";
        const bool sparse_shardmerge_mode = mode == "LEDGER_SPARSE_SHARDMERGE";
        const bool range_merge_mode = mode == "LEDGER_SHARDMERGE" ||
            persistent_shardmerge_mode || sparse_shardmerge_mode || hot_cache_mode;
        const bool atomic_shared_mode = mode == "ATOMIC_SHARED";
        const bool owner_tile_mode = mode == "OWNER_TILE";
        const bool batch_owner_mode = mode == "BATCH_OWNER";
        const bool cluster_pool_mode = mode == "LEDGER_CLUSTER_POOL";
        const bool adaptive_hot_cache_mode = mode == "LEDGER_ADAPTIVE_HOTCACHE";
        const bool adaptive_pool_mode = mode == "LEDGER_ADAPTIVE_POOL" || adaptive_hot_cache_mode;
        const bool auto_ledger_mode = mode == "LEDGER_AUTO";
        const bool cost_aware_mode = mode == "COST_AWARE" || mode == "LEDGER_AUTO_COST";
        const bool auto_cost_aware_mode = mode == "LEDGER_AUTO_COST";
        const bool chunked_mode = mode.rfind("CHUNKED_", 0) == 0;
        const bool block_pool_mode = mode == "LEDGER_BLOCK_POOL" ||
            mode.rfind("LEDGER_BLOCK_POOL_", 0) == 0;
        size_t block_pool_threshold = 256U;
        if (mode.rfind("LEDGER_BLOCK_POOL_", 0) == 0) {
            block_pool_threshold = parse_positive_size(
                std::string_view(mode).substr(std::string_view("LEDGER_BLOCK_POOL_").size()),
                "block-pool threshold");
        }
        const bool ledger_mode = mode == "LEDGER_REDUCE" || mode.rfind("LEDGER_REDUCE_", 0) == 0 || auto_ledger_mode || auto_cost_aware_mode;
        size_t ledger_block_size = 128U;
        size_t dispatch_chunk_size = 1U;
        if (chunked_mode) {
            dispatch_chunk_size = parse_positive_size(
                std::string_view(mode).substr(std::string_view("CHUNKED_").size()),
                "dispatch chunk size");
        }
        if (mode.rfind("LEDGER_REDUCE_", 0) == 0) {
            ledger_block_size = parse_positive_size(
                std::string_view(mode).substr(std::string_view("LEDGER_REDUCE_").size()),
                "ledger block size");
        }
        if (mode.rfind("SHARDED_", 0) == 0) {
            ledger_block_size = parse_positive_size(
                std::string_view(mode).substr(std::string_view("SHARDED_").size()),
                "sharded block size");
        }
        if (mode != "NO_STOP" && !ledger_mode && !sharded_mode &&
            !range_merge_mode && !atomic_shared_mode && !owner_tile_mode &&
            !batch_owner_mode && !cluster_pool_mode && !adaptive_pool_mode &&
            !auto_ledger_mode && !cost_aware_mode && !chunked_mode && !block_pool_mode) {
            throw std::runtime_error("unsupported mode");
        }
        if (mu_star != 0.0 || delta_flip != 0.0) throw std::runtime_error("NO_STOP requires zero stop parameters");
        if (!fs::is_directory(output_dir)) throw std::runtime_error("output directory missing");
        const fs::path cpu_tmp = output_dir / ".CPU_TOPOLOGY_AUDIT.tsv.tmp";
        const fs::path pred_tmp = output_dir / ".pred.tsv.tmp";
        const fs::path stops_tmp = output_dir / ".STOP_POINT_WALL.tsv.tmp";
        const fs::path workers_tmp = output_dir / ".WORKER_LEDGER.tsv.tmp";
        const fs::path metrics_tmp = output_dir / ".RUNNER_METRICS.tsv.tmp";
        const fs::path marker_tmp = output_dir / ".RUNNER_COMPLETE.json.tmp";
        for (const char* name : {
                 "CPU_TOPOLOGY_AUDIT.tsv", "pred.tsv", "STOP_POINT_WALL.tsv",
                 "WORKER_LEDGER.tsv", "RUNNER_METRICS.tsv", "RUNNER_COMPLETE.json",
                 ".CPU_TOPOLOGY_AUDIT.tsv.tmp", ".pred.tsv.tmp", ".STOP_POINT_WALL.tsv.tmp",
                 ".WORKER_LEDGER.tsv.tmp", ".RUNNER_METRICS.tsv.tmp", ".RUNNER_COMPLETE.json.tmp"}) {
            if (fs::exists(output_dir / name)) throw std::runtime_error("refuse existing output artifact");
        }

        pin_current_thread(control_cpu);
        const bool allow_hyperthreads = workers > 20U;
        write_cpu_audit(cpu_tmp, control_cpu, worker_cpus, allow_hyperthreads);
        acor_edit_distance_self_test();

        rusage usage_before{};
        getrusage(RUSAGE_SELF, &usage_before);
        const std::string load_t0 = load_average();
        const auto total_started = std::chrono::steady_clock::now();
        size_t total_reads = 0;
        std::vector<ClusterTask> tasks = load_tasks(
            reads_path, order_path, dataset, config_k(config_name), total_reads);
        const auto input_loaded = std::chrono::steady_clock::now();

        std::vector<ClusterResult> results(tasks.size());
        std::vector<WorkerStats> worker_stats(workers);
        const auto compute_started = std::chrono::steady_clock::now();
        double worker_start_seconds = 0.0;
        auto run_cluster_worker_pool = [&](const std::vector<size_t>& indices) {
            const auto pool_started = std::chrono::steady_clock::now();
            std::atomic<size_t> next{0U};
            std::atomic<size_t> ready{0U};
            std::atomic<bool> start{false};
            std::atomic<bool> failed{false};
            std::exception_ptr worker_error;
            std::mutex error_mutex;
            auto worker_loop = [&](size_t worker) {
                bool announced_ready = worker == 0U;
                try {
                    pin_current_thread(worker_cpus[worker]);
                    // Allocate and first-touch the worker's dense evidence
                    // ledger only after pinning this thread.  Constructing all
                    // ledgers on the control CPU places their pages on one
                    // NUMA node and makes workers on the other socket pay for
                    // remote memory throughout accumulation and decoding.
                    WorkerContext context;
                    if (worker != 0U) {
                        ready.fetch_add(1U, std::memory_order_release);
                        announced_ready = true;
                        while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                    }
                    while (!failed.load(std::memory_order_relaxed)) {
                        const size_t begin = next.fetch_add(dispatch_chunk_size, std::memory_order_relaxed);
                        if (begin >= indices.size()) break;
                        const size_t end = std::min(indices.size(), begin + dispatch_chunk_size);
                        for (size_t position = begin; position < end; ++position) {
                            const size_t index = indices[position];
                            results[index] = process_cluster(
                                tasks[index], context, worker, profile_enabled,
                                adaptive_hot_cache_mode);
                            ++worker_stats[worker].clusters;
                            worker_stats[worker].reads_consumed += results[index].stop_index;
                            worker_stats[worker].cluster_seconds_sum += results[index].cluster_seconds;
                        }
                    }
                } catch (...) {
                    failed.store(true, std::memory_order_relaxed);
                    if (!announced_ready) ready.fetch_add(1U, std::memory_order_release);
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!worker_error) worker_error = std::current_exception();
                }
            };
            std::vector<std::thread> threads;
            threads.reserve(workers > 0U ? workers - 1U : 0U);
            ThreadJoiner thread_joiner(threads);
            for (size_t worker = 1; worker < workers; ++worker) {
                threads.emplace_back([&, worker]() { worker_loop(worker); });
            }
            while (ready.load(std::memory_order_acquire) + 1U != workers) std::this_thread::yield();
            const auto workers_ready = std::chrono::steady_clock::now();
            worker_start_seconds += std::chrono::duration<double>(
                workers_ready - pool_started).count();
            start.store(true, std::memory_order_release);
            worker_loop(0U);
            for (std::thread& thread : threads) thread.join();
            if (worker_error) std::rethrow_exception(worker_error);
        };
        if (cost_aware_mode && !auto_cost_aware_mode) {
            // Cost-aware scheduling keeps the original cluster-level
            // reconstruction path but dispatches the largest independent
            // clusters first.  This is deliberately an opt-in mode: it is a
            // scheduling experiment, not an algorithmic change.
            std::vector<size_t> cost_order(tasks.size());
            std::iota(cost_order.begin(), cost_order.end(), 0U);
            sort_indices_by_estimated_work(cost_order, tasks);
            run_cluster_worker_pool(cost_order);
        } else if (cluster_pool_mode || adaptive_pool_mode) {
            // Persistent cluster-level pool for ledger reconstruction. Each
            // worker owns one CPU and reuses the scheduling loop across all
            // clusters; a cluster itself is processed with one dense ledger
            // worker, so this mode avoids nested oversubscription while
            // preserving the exact evidence path and deterministic output.
            std::vector<size_t> pooled_indices;
            pooled_indices.reserve(tasks.size());
            if (adaptive_pool_mode) {
                // Very large clusters benefit from the existing intra-cluster
                // range merge; sending them through a one-worker pool would
                // throw away useful parallelism.  The threshold is deliberately
                // conservative and is recorded by the benchmark mode.
                const size_t large_cluster_threshold = std::max<size_t>(10000U, workers * 512U);
                for (size_t index = 0; index < tasks.size(); ++index) {
                    if (tasks[index].reads.size() >= large_cluster_threshold) {
                        results[index] = process_ledger_cluster(
                            tasks[index], workers, worker_cpus, profile_enabled,
                            worker_stats, true, adaptive_hot_cache_mode);
                    } else {
                        pooled_indices.push_back(index);
                    }
                }
            } else {
                pooled_indices.resize(tasks.size());
                std::iota(pooled_indices.begin(), pooled_indices.end(), 0U);
            }
            // Small and medium clusters belong on the persistent cluster
            // pool.  The previous prototype called process_ledger_cluster()
            // with one worker for every such cluster, which recreated dense
            // ledgers and performed a one-way merge thousands of times.  That
            // was exact but defeated both state reuse and cache locality.
            // Reuse one model per pinned worker instead; process_cluster()
            // follows the same deterministic final-only evidence/decode path.
            run_cluster_worker_pool(pooled_indices);
        } else if (block_pool_mode) {
            // Persistent inner-ledger pool.  Small clusters stay on the
            // ordinary cluster-level path; only clusters above the explicit
            // threshold are split into deterministic read blocks.  The pool
            // is created once for all eligible clusters, avoiding one thread
            // creation/synchronization cycle per large cluster.
            std::vector<size_t> small_indices;
            std::vector<size_t> large_indices;
            small_indices.reserve(tasks.size());
            large_indices.reserve(tasks.size());
            // Splitting every cluster above the nominal block size is too
            // aggressive for small-read workloads: it creates many local
            // ledgers without exposing useful inner parallelism.  Restrict
            // the experimental path to genuinely large clusters and keep a
            // single-worker run identical to the baseline.
            const size_t large_cluster_threshold = std::max<size_t>(
                4096U, block_pool_threshold * 4U);
            for (size_t index = 0; index < tasks.size(); ++index) {
                if (workers > 1U && tasks[index].reads.size() >= large_cluster_threshold) {
                    large_indices.push_back(index);
                } else {
                    small_indices.push_back(index);
                }
            }
            run_cluster_worker_pool(small_indices);
            run_persistent_ledger_block_pool(
                large_indices, tasks, results, workers, worker_cpus,
                profile_enabled, block_pool_threshold, worker_stats);
        } else if (persistent_shardmerge_mode) {
            // Adaptive persistent shard-merge policy.  Medium and small
            // clusters stay on the persistent cluster-level pool, while only
            // genuinely large clusters pay for intra-cluster ledgers.  This
            // keeps the exact same evidence path but avoids multiplying dense
            // ledgers and merge work for thousands of modest clusters.
            const size_t large_cluster_threshold = std::max<size_t>(
                4096U, workers * 1024U);
            std::vector<size_t> small_indices;
            std::vector<size_t> large_indices;
            small_indices.reserve(tasks.size());
            large_indices.reserve(tasks.size());
            for (size_t index = 0; index < tasks.size(); ++index) {
                if (workers > 1U && tasks[index].reads.size() >= large_cluster_threshold) {
                    large_indices.push_back(index);
                } else {
                    small_indices.push_back(index);
                }
            }
            run_cluster_worker_pool(small_indices);
            for (const size_t index : large_indices) {
                results[index] = process_ledger_cluster(
                    tasks[index], workers, worker_cpus, profile_enabled,
                    worker_stats, true, false, true);
            }
        } else if (auto_ledger_mode || auto_cost_aware_mode) {
            // Keep ordinary clusters on the persistent cluster-level pool and
            // invoke intra-cluster ledger parallelism only for genuinely large
            // clusters. This prevents LEDGER_AUTO from serialising small jobs.
            const size_t ledger_threshold = std::max<size_t>(2048U, workers * 1024U);
            std::vector<size_t> small_indices;
            std::vector<size_t> large_indices;
            small_indices.reserve(tasks.size());
            for (size_t index = 0; index < tasks.size(); ++index) {
                if (tasks[index].reads.size() >= ledger_threshold) {
                    large_indices.push_back(index);
                } else {
                    small_indices.push_back(index);
                }
            }
            if (auto_cost_aware_mode) {
                sort_indices_by_estimated_work(small_indices, tasks);
                sort_indices_by_estimated_work(large_indices, tasks);
            }
            run_cluster_worker_pool(small_indices);
            for (const size_t index : large_indices) {
                results[index] = process_ledger_cluster(
                    tasks[index], workers, worker_cpus, profile_enabled,
                    worker_stats, true, false);
            }
        } else if (ledger_mode || range_merge_mode || atomic_shared_mode || owner_tile_mode ||
                   batch_owner_mode) {
            // Hybrid policy: retain the existing cluster-level path for small
            // clusters, and expose intra-cluster parallelism only when a
            // cluster contains at least two requested ledger blocks.  This
            // avoids paying dense-ledger construction/merge costs on tiny
            // tasks while removing the dominant heavy-tail bottleneck.
            WorkerContext fallback_context;
            const size_t ledger_threshold = auto_ledger_mode ? std::max<size_t>(2048U, workers * 1024U) : ledger_block_size * 2U;
            for (size_t cluster_index = 0; cluster_index < tasks.size(); ++cluster_index) {
                const ClusterTask& task = tasks[cluster_index];
                if (task.reads.size() >= ledger_threshold) {
                    if (atomic_shared_mode) {
                        results[cluster_index] = process_atomic_shared_cluster(
                            task, workers, worker_cpus, profile_enabled, worker_stats);
                    } else if (owner_tile_mode) {
                        results[cluster_index] = process_owner_tile_cluster(
                            task, workers, worker_cpus, profile_enabled, worker_stats);
                    } else if (batch_owner_mode) {
                        results[cluster_index] = process_batch_owner_cluster(
                            task, workers, worker_cpus, profile_enabled, worker_stats);
                    } else {
                        results[cluster_index] = sharded_mode
                            ? process_sharded_cluster(
                                task, workers, worker_cpus, profile_enabled, worker_stats)
                            : process_ledger_cluster(
                                task, workers, worker_cpus, profile_enabled, worker_stats,
                                range_merge_mode || auto_ledger_mode, hot_cache_mode,
                                persistent_shardmerge_mode, sparse_shardmerge_mode);
                    }
                } else {
                    results[cluster_index] = process_cluster(
                        task, fallback_context, 0U, profile_enabled);
                    ++worker_stats[0].clusters;
                    worker_stats[0].reads_consumed += results[cluster_index].stop_index;
                    worker_stats[0].cluster_seconds_sum += results[cluster_index].cluster_seconds;
                }
            }
        } else {
            std::vector<size_t> all_indices(tasks.size());
            std::iota(all_indices.begin(), all_indices.end(), 0U);
            if (cost_aware_mode) sort_indices_by_estimated_work(all_indices, tasks);
            run_cluster_worker_pool(all_indices);
        }
        pin_current_thread(control_cpu);
        const auto compute_finished = std::chrono::steady_clock::now();

        size_t consumed_reads = 0;
        size_t early_clusters = 0;
        size_t completed_empty = 0;
        double graph_seconds = 0.0;
        double beam_seconds = 0.0;
        double evidence_seconds = 0.0;
        double model_reset_seconds = 0.0;
        double canonical_seconds = 0.0;
        double add_read_seconds = 0.0;
        double soft_decode_seconds = 0.0;
        double cluster_seconds_sum = 0.0;
        double ledger_accumulation_seconds = 0.0;
        double ledger_merge_seconds = 0.0;
        double ledger_final_decode_seconds = 0.0;
        std::ofstream pred(pred_tmp);
        pred << "row_key\tdataset_id\tcluster_id\tlegacy_prefix_alias\tstop_index"
                "\tconsumed_reads\ttotal_reads\treconstructed_sequence\tstatus\n";
        for (const ClusterResult& result : results) {
            consumed_reads += result.stop_index;
            early_clusters += static_cast<size_t>(result.stopped_early);
            completed_empty += static_cast<size_t>(result.completed_empty);
            graph_seconds += result.graph_update_seconds;
            beam_seconds += result.beam_decode_seconds;
            evidence_seconds += result.evidence_ed_seconds;
            cluster_seconds_sum += result.cluster_seconds;
            model_reset_seconds += result.model_reset_seconds;
            canonical_seconds += result.canonical_seconds;
            add_read_seconds += result.add_read_seconds;
            soft_decode_seconds += result.soft_decode_seconds;
            ledger_accumulation_seconds += result.ledger_accumulation_seconds;
            ledger_merge_seconds += result.ledger_merge_seconds;
            ledger_final_decode_seconds += result.ledger_final_decode_seconds;
            pred << output_dataset_id << '|' << result.cluster_id << '|' << result.legacy_prefix_alias
                 << '\t' << output_dataset_id << '\t' << result.cluster_id << '\t'
                 << result.legacy_prefix_alias << '\t' << result.stop_index << '\t'
                 << result.stop_index << '\t' << result.total_reads << '\t'
                 << result.candidate << '\t'
                 << (result.completed_empty ? "completed_empty" : "completed_nonempty") << '\n';
        }
        pred.flush();
        if (!pred) throw std::runtime_error("pred write failed");
        const auto timed_output_finished = std::chrono::steady_clock::now();
        pred.close();
        if (!pred) throw std::runtime_error("pred close failed");
        const std::string load_t1 = load_average();

        std::ofstream stops(stops_tmp);
        stops << "dataset\tround_id\tcluster_id\tstop_index\ttotal_reads\tunused_reads\tstopped_early\toutput_status\n";
        for (const ClusterResult& result : results) {
            stops << dataset << '\t' << round_id << '\t' << result.cluster_id << '\t'
                  << result.stop_index << '\t' << result.total_reads << '\t'
                  << (result.total_reads - result.stop_index) << '\t'
                  << (result.stopped_early ? 1 : 0) << '\t'
                  << (result.completed_empty ? "completed_empty" : "completed_nonempty") << '\n';
        }
        stops.close();
        if (!stops) throw std::runtime_error("stop write failed");

        std::ofstream worker_output(workers_tmp);
        worker_output << "worker_id\tlogical_cpu\tclusters\treads_consumed\tcluster_seconds_sum\n";
        worker_output << std::setprecision(12);
        for (size_t worker = 0; worker < workers; ++worker) {
            worker_output << worker << '\t' << worker_cpus[worker] << '\t'
                          << worker_stats[worker].clusters << '\t'
                          << worker_stats[worker].reads_consumed << '\t'
                          << worker_stats[worker].cluster_seconds_sum << '\n';
        }
        worker_output.flush();
        worker_output.close();
        if (!worker_output) throw std::runtime_error("worker ledger write/close failed");

        rusage usage_after{};
        getrusage(RUSAGE_SELF, &usage_after);
        const double input_seconds = std::chrono::duration<double>(input_loaded - total_started).count();
        const double thread_start_seconds = worker_start_seconds;
        const double compute_seconds = std::chrono::duration<double>(compute_finished - compute_started).count();
        const double parallel_active_seconds = std::max(
            0.0, compute_seconds - thread_start_seconds);
        const double pred_write_seconds = std::chrono::duration<double>(timed_output_finished - compute_finished).count();
        const double total_seconds = std::chrono::duration<double>(timed_output_finished - total_started).count();
        std::set<std::pair<int, int>> physical_worker_cores;
        for (int cpu : worker_cpus) {
            const CpuLocation location = cpu_location(cpu);
            physical_worker_cores.insert({location.socket, location.core});
        }
        std::ofstream metrics(metrics_tmp);
        metrics << "dataset\tround_id\tworkers\tworker_cpus\tcontrol_cpu\tphysical_cores\thyperthread_diagnostic\tprofile_enabled\tclusters\ttotal_reads\tconsumed_reads"
                   "\treads_saved_fraction\tearly_clusters\tcompleted_empty\tinput_load_seconds"
                   "\tthread_start_seconds\tcompute_seconds\tparallel_active_seconds\tpred_write_seconds\ttotal_wall_seconds"
                   "\tmodel_reset_cpu_seconds_sum\tcanonical_cpu_seconds_sum\tadd_read_cpu_seconds_sum"
                   "\tsoft_decode_cpu_seconds_sum\tgraph_update_cpu_seconds_sum\tbeam_decode_cpu_seconds_sum\tevidence_ed_cpu_seconds_sum"
                   "\tcluster_cpu_seconds_sum\tledger_accumulation_wall_seconds\tledger_merge_wall_seconds\tledger_final_decode_wall_seconds"
                   "\tuser_cpu_seconds\tsystem_cpu_seconds\tmax_rss_kb"
                   "\tloadavg_t0\tloadavg_t1\tfuture_read_violations\tsafe_batch_violations"
                   "\tstop_index_violations\tstatus\n";
        const bool run_pass = completed_empty == 0U && consumed_reads == total_reads;
        metrics << std::setprecision(15)
                << dataset << '\t' << round_id << '\t' << workers << '\t' << argv[8]
                << '\t' << control_cpu << '\t'
                << physical_worker_cores.size()
                << '\t' << (allow_hyperthreads ? 1 : 0) << '\t' << (profile_enabled ? 1 : 0)
                << '\t' << results.size() << '\t' << total_reads << '\t'
                << consumed_reads << "\tNOT_APPLICABLE\tNOT_APPLICABLE\t"
                << completed_empty << '\t'
                << input_seconds << '\t' << thread_start_seconds << '\t' << compute_seconds << '\t'
                << parallel_active_seconds << '\t' << pred_write_seconds << '\t' << total_seconds << '\t'
                << model_reset_seconds << '\t'
                << canonical_seconds << '\t' << add_read_seconds << '\t' << soft_decode_seconds << '\t'
                << graph_seconds << '\t'
                << beam_seconds << '\t' << evidence_seconds << '\t' << cluster_seconds_sum << '\t'
                << ledger_accumulation_seconds << '\t' << ledger_merge_seconds << '\t'
                << ledger_final_decode_seconds << '\t'
                << (timeval_seconds(usage_after.ru_utime) - timeval_seconds(usage_before.ru_utime)) << '\t'
                << (timeval_seconds(usage_after.ru_stime) - timeval_seconds(usage_before.ru_stime)) << '\t'
                << usage_after.ru_maxrss << '\t' << load_t0 << '\t' << load_t1
                << "\tNOT_APPLICABLE\tNOT_APPLICABLE\tNOT_APPLICABLE\t"
                << (run_pass ? "PASS" : "FAIL") << '\n';
        metrics.close();
        if (!metrics) throw std::runtime_error("metrics write failed");
        fs::rename(cpu_tmp, output_dir / "CPU_TOPOLOGY_AUDIT.tsv");
        fs::rename(pred_tmp, output_dir / "pred.tsv");
        fs::rename(stops_tmp, output_dir / "STOP_POINT_WALL.tsv");
        fs::rename(workers_tmp, output_dir / "WORKER_LEDGER.tsv");
        fs::rename(metrics_tmp, output_dir / "RUNNER_METRICS.tsv");
        {
            std::ofstream marker(marker_tmp);
            marker << "{\"state\":\"COMPLETE\",\"status\":\""
                   << (run_pass ? "PASS" : "FAIL") << "\"}\n";
            marker.flush();
            if (!marker) throw std::runtime_error("runner completion marker write failed");
            marker.close();
            if (!marker) throw std::runtime_error("runner completion marker close failed");
        }
        fs::rename(marker_tmp, output_dir / "RUNNER_COMPLETE.json");
        std::cout << std::setprecision(15) << total_seconds << '\n';
        return run_pass ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
