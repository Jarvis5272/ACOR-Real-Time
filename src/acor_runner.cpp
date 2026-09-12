#include "acor_kernel.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <pthread.h>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <sys/resource.h>
#include <thread>
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
};

struct WorkerStats {
    size_t clusters = 0;
    size_t reads_consumed = 0;
    double cluster_seconds_sum = 0.0;
};

struct WorkerContext {
    AcorKernel model;
};

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
    double value = 0.0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc() || parsed.ptr != text.data() + text.size() || !std::isfinite(value)) {
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

static ClusterResult process_cluster(
    const ClusterTask& task, WorkerContext& context, size_t worker_id,
    bool profile_enabled) {
    const auto cluster_started = std::chrono::steady_clock::now();
    auto profile_started = cluster_started;
    double model_reset_seconds = 0.0;
    if (profile_enabled) profile_started = std::chrono::steady_clock::now();
    context.model.reset(task.design_length);
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

    for (size_t index = 0; index < task.reads.size(); ++index) {
        const bool is_last = index + 1U == task.reads.size();
        if (profile_enabled) profile_started = std::chrono::steady_clock::now();
        const std::string sequence = acor_canonical(task.reads[index]);
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
        if (config_name != "kmc_k9_b16_lognormal") throw std::runtime_error("unsupported frozen config");
        if (mode != "NO_STOP") throw std::runtime_error("only NO_STOP is supported");
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
            reads_path, order_path, dataset, 9, total_reads);
        const auto input_loaded = std::chrono::steady_clock::now();

        std::vector<ClusterResult> results(tasks.size());
        std::vector<WorkerStats> worker_stats(workers);
        std::atomic<size_t> next{0U};
        std::atomic<size_t> ready{0U};
        std::atomic<bool> start{false};
        std::atomic<bool> failed{false};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        std::vector<std::unique_ptr<WorkerContext>> contexts;
        contexts.reserve(workers);
        for (size_t worker = 0; worker < workers; ++worker) {
            contexts.push_back(std::make_unique<WorkerContext>());
        }
        auto worker_loop = [&](size_t worker) {
            bool announced_ready = worker == 0U;
            try {
                pin_current_thread(worker_cpus[worker]);
                if (worker != 0U) {
                    ready.fetch_add(1U, std::memory_order_release);
                    announced_ready = true;
                    while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                }
                while (!failed.load(std::memory_order_relaxed)) {
                    const size_t index = next.fetch_add(1U, std::memory_order_relaxed);
                    if (index >= tasks.size()) break;
                    results[index] = process_cluster(
                        tasks[index], *contexts[worker], worker, profile_enabled);
                    ++worker_stats[worker].clusters;
                    worker_stats[worker].reads_consumed += results[index].stop_index;
                    worker_stats[worker].cluster_seconds_sum += results[index].cluster_seconds;
                }
            } catch (...) {
                failed.store(true, std::memory_order_relaxed);
                if (!announced_ready) {
                    ready.fetch_add(1U, std::memory_order_release);
                }
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };
        std::vector<std::thread> threads;
        threads.reserve(workers > 0U ? workers - 1U : 0U);
        ThreadJoiner thread_joiner(threads);
        for (size_t worker = 1; worker < workers; ++worker) {
            threads.emplace_back([&, worker]() {
                worker_loop(worker);
            });
        }
        while (ready.load(std::memory_order_acquire) + 1U != workers) std::this_thread::yield();
        const auto compute_started = std::chrono::steady_clock::now();
        start.store(true, std::memory_order_release);
        worker_loop(0U);
        for (std::thread& thread : threads) thread.join();
        pin_current_thread(control_cpu);
        const auto compute_finished = std::chrono::steady_clock::now();
        if (worker_error) std::rethrow_exception(worker_error);

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
        const double thread_start_seconds = std::chrono::duration<double>(compute_started - input_loaded).count();
        const double compute_seconds = std::chrono::duration<double>(compute_finished - compute_started).count();
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
                   "\tthread_start_seconds\tcompute_seconds\tpred_write_seconds\ttotal_wall_seconds"
                   "\tmodel_reset_cpu_seconds_sum\tcanonical_cpu_seconds_sum\tadd_read_cpu_seconds_sum"
                   "\tsoft_decode_cpu_seconds_sum\tgraph_update_cpu_seconds_sum\tbeam_decode_cpu_seconds_sum\tevidence_ed_cpu_seconds_sum"
                   "\tcluster_cpu_seconds_sum\tuser_cpu_seconds\tsystem_cpu_seconds\tmax_rss_kb"
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
                << pred_write_seconds << '\t' << total_seconds << '\t' << model_reset_seconds << '\t'
                << canonical_seconds << '\t' << add_read_seconds << '\t' << soft_decode_seconds << '\t'
                << graph_seconds << '\t'
                << beam_seconds << '\t' << evidence_seconds << '\t' << cluster_seconds_sum << '\t'
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
