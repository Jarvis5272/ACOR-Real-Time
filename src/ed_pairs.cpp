#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

static size_t edit_distance(const std::string& lhs, const std::string& rhs) {
    const std::string& rows = lhs.size() >= rhs.size() ? lhs : rhs;
    const std::string& columns = lhs.size() >= rhs.size() ? rhs : lhs;
    std::vector<size_t> previous(columns.size() + 1U);
    std::vector<size_t> current(columns.size() + 1U);
    for (size_t j = 0; j <= columns.size(); ++j) previous[j] = j;
    for (size_t i = 1; i <= rows.size(); ++i) {
        current[0] = i;
        for (size_t j = 1; j <= columns.size(); ++j) {
            const size_t cost = rows[i - 1U] == columns[j - 1U] ? 0U : 1U;
            current[j] = std::min({
                previous[j] + 1U,
                current[j - 1U] + 1U,
                previous[j - 1U] + cost,
            });
        }
        previous.swap(current);
    }
    return previous.back();
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        if (line.back() == '\r') {
            std::cerr << "CRLF input is not accepted\n";
            return 2;
        }
        const size_t first = line.find('\t');
        if (first == std::string::npos) {
            std::cerr << "missing key delimiter: " << line << '\n';
            return 2;
        }
        const size_t second = line.find('\t', first + 1);
        if (second == std::string::npos) {
            std::cerr << "missing second delimiter: " << line << '\n';
            return 2;
        }
        if (line.find('\t', second + 1) != std::string::npos) {
            std::cerr << "too many columns: " << line << '\n';
            return 2;
        }
        const std::string key = line.substr(0, first);
        if (key.empty()) {
            std::cerr << "empty key\n";
            return 2;
        }
        const std::string a = line.substr(first + 1, second - first - 1);
        const std::string b = line.substr(second + 1);
        std::cout << key << '\t' << a.size() << '\t' << b.size() << '\t'
                  << edit_distance(a, b) << '\n';
    }
    if (!std::cin.eof()) {
        std::cerr << "input stream ended without clean EOF\n";
        return 2;
    }
    std::cout.flush();
    if (!std::cout) {
        std::cerr << "output write failed\n";
        return 2;
    }
    return 0;
}
