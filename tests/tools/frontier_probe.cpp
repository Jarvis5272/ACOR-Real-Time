#define ACOR_FRONTIER_AUDIT
#define main acor_kernel_standalone_main_unused
#include "../../src/acor_kernel.cpp"
#undef main

#include <fstream>
#include <map>

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    std::ifstream input(argv[1]);
    std::ofstream frontier(argv[2]);
    std::ofstream final_output(argv[3]);
    if (!input || !frontier || !final_output) return 2;
    std::string line;
    if (!std::getline(input, line) || line != "cluster_id\tdesign_length\tread_sequence") return 2;
    struct Group { int design_length = 0; std::vector<std::string> reads; };
    std::map<std::string, Group> groups;
    std::vector<std::string> order;
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = split_tsv(line);
        if (fields.size() != 3U) return 2;
        if (groups.find(fields[0]) == groups.end()) order.push_back(fields[0]);
        Group& group = groups[fields[0]];
        const int design_length = std::stoi(fields[1]);
        if (group.design_length == 0) group.design_length = design_length;
        if (group.design_length != design_length) return 2;
        group.reads.push_back(fields[2]);
    }
    set_frontier_audit_output(&frontier);
    final_output << "cluster_id\treconstructed_sequence\n";
    const Config config = parse_config("kmc_k9_b16_lognormal");
    for (const std::string& cluster_id : order) {
        IncrementalKmcBeam model(config, ConstraintMode::SoftLength, OptimizationMode::FinalOnlyExact, 0);
        const Group& group = groups.at(cluster_id);
        model.reset(group.design_length);
        ConsensusResult result;
        for (size_t index = 0; index < group.reads.size(); ++index) {
            result = model.add_read(group.reads[index], index + 1U == group.reads.size());
        }
        final_output << cluster_id << '\t'
                     << (result.ranked.empty() ? std::string() : result.ranked.front().sequence)
                     << '\n';
    }
    return frontier && final_output ? 0 : 2;
}
