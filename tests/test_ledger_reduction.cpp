#include "acor_kernel.hpp"

#include <iostream>
#include <string>
#include <vector>

static std::string make_read(size_t seed, size_t length) {
    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    std::string sequence(length, 'A');
    for (size_t index = 0; index < length; ++index) {
        sequence[index] = bases[(index * 7U + seed * 11U + (index / 9U)) % 4U];
    }
    if (seed % 5U == 0U && length > 17U) sequence[17] = 'A';
    if (seed % 7U == 0U && length > 43U) sequence[43] = 'N';
    return sequence;
}

static std::string top_sequence(const AcorIncrementalResult& result) {
    return result.ranked.empty() ? std::string() : result.ranked.front().sequence;
}

int main() {
    constexpr int design_length = 100;
    constexpr size_t read_count = 96;
    std::vector<std::string> reads;
    reads.reserve(read_count);
    for (size_t index = 0; index < read_count; ++index) {
        reads.push_back(make_read(index, design_length));
    }

    AcorKernel serial;
    serial.reset(design_length);
    for (size_t index = 0; index < reads.size(); ++index) {
        (void)serial.add_read(reads[index], index + 1U == reads.size());
    }
    const std::string serial_top = top_sequence(serial.decode_final());

    AcorKernel merged;
    merged.reset(design_length);
    constexpr size_t blocks = 8;
    for (size_t block = 0; block < blocks; ++block) {
        AcorKernel local;
        local.reset(design_length);
        for (size_t index = block; index < reads.size(); index += blocks) {
            local.add_read_evidence_only(reads[index]);
        }
        merged.merge_evidence_from(local);
    }
    const std::string merged_top = top_sequence(merged.decode_final());

    if (serial_top != merged_top) {
        std::cerr << "ledger reduction mismatch\n"
                  << "serial=" << serial_top << "\n"
                  << "merged=" << merged_top << "\n";
        return 1;
    }
    std::cout << "PASS\tserial_top=" << serial_top
              << "\tmerged_top=" << merged_top << '\n';
    return 0;
}
