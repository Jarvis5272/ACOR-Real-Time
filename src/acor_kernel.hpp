#pragma once

#include <cstdint>
#include <memory>
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

class AcorKernel {
public:
    AcorKernel();
    ~AcorKernel();
    AcorKernel(AcorKernel&&) noexcept;
    AcorKernel& operator=(AcorKernel&&) noexcept;
    AcorKernel(const AcorKernel&) = delete;
    AcorKernel& operator=(const AcorKernel&) = delete;

    void reset(int design_length);
    AcorIncrementalResult add_read(const std::string& sequence, bool is_last);
    double graph_update_seconds() const;
    double beam_decode_seconds() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::string acor_canonical(const std::string& sequence);
void acor_edit_distance_self_test();
