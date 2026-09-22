#pragma once

#include <cstdint>
#include <memory>
#include <utility>
#include <string>
#include <vector>

struct AcorCandidate {
    std::string sequence;
    double score = 0.0;
};

struct AcorIncrementalResult {
    std::vector<AcorCandidate> ranked;
    double final_decode_seconds = 0.0;
};

// A counter update used by the owner-sharded evidence prototype.  The five
// counter kinds correspond to node, edge, outgoing, start, and end support.
// Updates are additive and therefore can be reduced in any fixed owner order.
struct AcorCounterDelta {
    std::uint8_t kind = 0;
    std::uint32_t index = 0;
    std::uint32_t amount = 0;
};

class AcorKernel {
public:
    AcorKernel();
    ~AcorKernel();
    AcorKernel(AcorKernel&&) noexcept;
    AcorKernel& operator=(AcorKernel&&) noexcept;
    AcorKernel(const AcorKernel&) = delete;
    AcorKernel& operator=(const AcorKernel&) = delete;

    void reset(int design_length);
    void reserve_evidence(std::size_t read_count);
    AcorIncrementalResult add_read(const std::string& sequence, bool is_last);
    void add_read_evidence_only(const std::string& sequence);
    void add_read_graph_atomic(const std::string& sequence);
    // Exact diagnostic path: cache repeated counter updates in a small
    // thread-local table, then flush before reduction/decode.
    void enable_hot_cache();
    void flush_pending_hot_cache();
    void enable_atomic_updates();
    void flush_atomic_updates();
    // Emit exactly the additive graph-counter updates produced by one read,
    // without touching the dense ledger.  This is used by the batched owner
    // reduction prototype; metadata is maintained by the caller.
    void append_evidence_deltas(
        const std::string& sequence, std::vector<AcorCounterDelta>& output) const;
    void merge_evidence_from(const AcorKernel& other);
    void merge_counter_range_from(
        const AcorKernel& other, std::size_t shard, std::size_t shard_count);
    // Sparse owner-sharded reduction: visit only non-zero counters owned by
    // this shard instead of scanning the full dense counter range.
    void merge_counter_sparse_range_from(
        const AcorKernel& other, std::size_t shard, std::size_t shard_count);
    void finalize_counter_range_merge();
    void merge_evidence_metadata_from(const AcorKernel& other);
    void apply_counter_deltas(const std::vector<AcorCounterDelta>& deltas);
    void apply_counter_tile(
        std::uint8_t kind, std::uint32_t begin, const std::uint32_t* counts,
        std::size_t count);
    void finalize_counter_deltas();
    void merge_sparse_metadata(
        const std::vector<std::pair<std::string, std::uint32_t>>& read_counts,
        const std::vector<int>& lengths,
        std::size_t reads_arrived,
        std::size_t exact_m_read_count,
        double log_length_sum,
        double log_length_sq_sum,
        std::uint64_t total_start_observations);
    AcorIncrementalResult decode_final();
    double graph_update_seconds() const;
    double beam_decode_seconds() const;
    bool uses_sparse_ledger() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::string acor_canonical(const std::string& sequence);
void acor_canonical_into(const std::string& sequence, std::string& output);
void acor_edit_distance_self_test();
