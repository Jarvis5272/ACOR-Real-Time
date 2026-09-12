#define main acor_kernel_standalone_main_unused
#include "../../src/acor_kernel.cpp"
#undef main

#include <iostream>

int main() {
    const double one = 1.0;
    const std::vector<Candidate> values = {
        {"same", 0.0}, {"same", -0.0}, {"same", 0.49e-12}, {"same", -0.49e-12},
        {"same", 0.51e-12}, {"same", -0.51e-12}, {"same", one},
        {"same", std::nextafter(one, 2.0)}, {"same", one + 0.49e-12},
        {"same", one + 0.51e-12}, {"alpha", one}, {"beta", one},
    };
    auto less = [](const Candidate& lhs, const Candidate& rhs) { return candidate_better(lhs, rhs); };
    auto equivalent = [&](const Candidate& lhs, const Candidate& rhs) {
        return !less(lhs, rhs) && !less(rhs, lhs);
    };
    for (const Candidate& a : values) {
        if (less(a, a)) return 1;
        for (const Candidate& b : values) {
            if (less(a, b) && less(b, a)) return 2;
            for (const Candidate& c : values) {
                if (less(a, b) && less(b, c) && !less(a, c)) return 3;
                if (equivalent(a, b) && equivalent(b, c) && !equivalent(a, c)) return 4;
            }
        }
    }
    int rejected = 0;
    for (double score : {std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()}) {
        try { (void) candidate_better(Candidate{"a", score}, Candidate{"b", 0.0}); }
        catch (const std::exception&) { ++rejected; }
    }
    return rejected == 3 ? 0 : 5;
}
