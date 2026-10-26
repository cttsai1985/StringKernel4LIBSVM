#include "StringKernels.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Feature = std::tuple<std::size_t, std::size_t, std::string>;
using FeatureMap = std::set<Feature>;

Catch::Approx approximately(double value) { return Catch::Approx(value).epsilon(1e-12); }

// This intentionally allocates explicit features, rather than following the
// production matching-prefix algorithm or calling its weighting helper.
FeatureMap features(std::string_view sequence, int degree) {
    FeatureMap result;
    const auto maximum = std::min(sequence.size(), static_cast<std::size_t>(degree));
    for (std::size_t order = 1; order <= maximum; ++order) {
        for (std::size_t position = 0; position <= sequence.size() - order; ++position) {
            result.emplace(position, order, std::string(sequence.substr(position, order)));
        }
    }
    return result;
}

double oracle_raw(std::string_view a, std::string_view b, int degree) {
    const auto left = features(a, degree);
    const auto right = features(b, degree);
    long double numerator = 0;
    for (const auto& feature : left) {
        if (right.count(feature) != 0) {
            numerator += static_cast<long double>(degree) - std::get<1>(feature) + 1;
        }
    }
    const auto d = static_cast<long double>(degree);
    return static_cast<double>(2 * numerator / (d * (d + 1)));
}

double oracle_normalized(std::string_view a, std::string_view b, int degree) {
    const auto left = oracle_raw(a, a, degree);
    const auto right = oracle_raw(b, b, degree);
    if (left == 0 || right == 0) {
        return 0;
    }
    return oracle_raw(a, b, degree) / std::sqrt(left * right);
}

std::vector<std::string> short_sequences() {
    std::vector<std::string> result{""};
    std::vector<std::string> level{""};
    for (int length = 1; length <= 3; ++length) {
        std::vector<std::string> next;
        for (const auto& prefix : level) {
            for (const auto symbol : std::string_view("AC")) {
                next.push_back(prefix + symbol);
            }
        }
        result.insert(result.end(), next.begin(), next.end());
        level = std::move(next);
    }
    return result;
}

using TokenCounts = std::unordered_map<std::string, long double>;

TokenCounts spectrum_features(std::string_view sequence, int k) {
    TokenCounts counts;
    if (k <= 0 || static_cast<std::size_t>(k) > sequence.size()) {
        return counts;
    }
    for (std::size_t i = 0; i + static_cast<std::size_t>(k) <= sequence.size(); ++i) {
        ++counts[std::string(sequence.substr(i, static_cast<std::size_t>(k)))];
    }
    return counts;
}

TokenCounts subsequence_features(std::string_view sequence, int k, double decay) {
    TokenCounts counts;
    if (k <= 0 || static_cast<std::size_t>(k) > sequence.size()) {
        return counts;
    }
    std::vector<std::size_t> positions(static_cast<std::size_t>(k));
    const auto visit = [&](const auto& self, std::size_t depth, std::size_t start) -> void {
        if (depth == positions.size()) {
            std::string token;
            for (const auto position : positions) {
                token.push_back(sequence[position]);
            }
            const auto gaps = positions.back() - positions.front() + 1 - positions.size();
            counts[token] += std::pow(static_cast<long double>(decay), gaps);
            return;
        }
        const auto remaining = positions.size() - depth;
        for (std::size_t position = start; position + remaining <= sequence.size(); ++position) {
            positions[depth] = position;
            self(self, depth + 1, position + 1);
        }
    };
    visit(visit, 0, 0);
    return counts;
}

double feature_dot(const TokenCounts& left, const TokenCounts& right) {
    long double result = 0;
    for (const auto& [token, count] : left) {
        const auto found = right.find(token);
        if (found != right.end()) {
            result += count * found->second;
        }
    }
    return static_cast<double>(result);
}

double oracle_feature_kernel(const TokenCounts& left, const TokenCounts& right) {
    const auto left_self = feature_dot(left, left);
    const auto right_self = feature_dot(right, right);
    if (left_self == 0 || right_self == 0) {
        return 0;
    }
    return feature_dot(left, right) / std::sqrt(left_self * right_self);
}

} // namespace

TEST_CASE("Default kernel settings remain weighted degree five", "[kernel][defaults]") {
    const stringkernel::KernelConfig config;
    CHECK(config.type == stringkernel::KernelType::weighted_degree);
    CHECK(config.degree == 5);
    CHECK(static_cast<int>(stringkernel::KernelType::distance) == 0);
    CHECK(static_cast<int>(stringkernel::KernelType::weighted_degree) == 1);
    CHECK(static_cast<int>(stringkernel::KernelType::spectrum) == 2);
    CHECK(static_cast<int>(stringkernel::KernelType::subsequence) == 3);
    CHECK(config.decay == approximately(0.5));
    CHECK_NOTHROW(config.validate());
}

TEST_CASE("Spectrum matches an independently counted fixed-k feature map", "[kernel][spectrum]") {
    using namespace stringkernel;
    for (const auto k : {1, 2, 3, 4}) {
        for (const auto& a : short_sequences()) {
            CAPTURE(k, a);
            const auto a_features = spectrum_features(a, k);
            CHECK(spectrum_kernel(a, a, k) ==
                  approximately(feature_dot(a_features, a_features)).margin(1e-13));
            CHECK(spectrum_self(a, k) ==
                  approximately(feature_dot(a_features, a_features)).margin(1e-13));
            for (const auto& b : short_sequences()) {
                CAPTURE(b);
                CHECK(
                    spectrum_kernel(a, b, k) ==
                    approximately(feature_dot(a_features, spectrum_features(b, k))).margin(1e-13));
                CHECK(kernel_value(a, b, {KernelType::spectrum, k}) ==
                      approximately(oracle_feature_kernel(a_features, spectrum_features(b, k)))
                          .margin(1e-13));
            }
        }
    }
    CHECK(spectrum_kernel("AAAA", "AAA", 2) == 6);
    CHECK(kernel_value("ACG", "GAC", {KernelType::spectrum, 2}) > 0);
    CHECK(kernel_value("ACG", "GAC", {KernelType::weighted_degree, 2}) == 0);
    CHECK(kernel_value("ACG", "AG", {KernelType::spectrum, 2}) == 0);
}

TEST_CASE("Subsequence kernel matches explicit ordered embedding features",
          "[kernel][subsequence]") {
    using namespace stringkernel;
    for (const auto k : {1, 2, 3, 4}) {
        for (const auto decay : {0.1, 0.5, 1.0}) {
            for (const auto& a : short_sequences()) {
                CAPTURE(k, decay, a);
                const auto a_features = subsequence_features(a, k, decay);
                CHECK(subsequence_kernel(a, a, k, decay) ==
                      approximately(feature_dot(a_features, a_features)).margin(1e-13));
                CHECK(subsequence_self(a, k, decay) ==
                      approximately(feature_dot(a_features, a_features)).margin(1e-13));
                for (const auto& b : short_sequences()) {
                    CAPTURE(b);
                    const auto b_features = subsequence_features(b, k, decay);
                    CHECK(subsequence_kernel(a, b, k, decay) ==
                          approximately(feature_dot(a_features, b_features)).margin(1e-13));
                    CHECK(
                        kernel_value(a, b, {KernelType::subsequence, k, decay}) ==
                        approximately(oracle_feature_kernel(a_features, b_features)).margin(1e-13));
                }
            }
        }
    }
    CHECK(kernel_value("ACG", "AG", {KernelType::subsequence, 2, 0.5}) == approximately(1.0 / 3.0));
    CHECK(subsequence_kernel("AAAA", "AAAA", 2, 1.0) == 36);
    CHECK(subsequence_self("AAAA", 5, 1.0) == 0);
}

TEST_CASE("Spectrum and subsequence reject invalid parameters and nonfinite arithmetic",
          "[kernel][validation]") {
    using namespace stringkernel;
    for (const auto k : {0, -1}) {
        CHECK_THROWS_AS(spectrum_kernel("A", "A", k), std::invalid_argument);
        CHECK_THROWS_AS(spectrum_self("A", k), std::invalid_argument);
        CHECK_THROWS_AS(subsequence_kernel("A", "A", k, 0.5), std::invalid_argument);
        CHECK_THROWS_AS(subsequence_self("A", k, 0.5), std::invalid_argument);
    }
    for (const auto decay : {0.0, -0.1, 1.1, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
        CHECK_THROWS_AS((KernelConfig{KernelType::subsequence, 2, decay}.validate()),
                        std::invalid_argument);
        CHECK_THROWS_AS(subsequence_kernel("ABC", "ABC", 2, decay), std::invalid_argument);
    }
    CHECK_THROWS_AS((KernelConfig{KernelType::spectrum, 2, 0.25}.validate()),
                    std::invalid_argument);
    CHECK_THROWS_AS((KernelConfig{KernelType::weighted_degree, 2, 0.25}.validate()),
                    std::invalid_argument);
    CHECK_THROWS_AS(subsequence_self(std::string(550, 'A'), 275, 1.0), std::overflow_error);
}

TEST_CASE("Tiny positive SSK decay retains nonempty gap-free diagonal features",
          "[kernel][subsequence][numeric]") {
    using namespace stringkernel;
    constexpr double tiny_decay = 1e-300;
    CHECK(subsequence_self("AA", 2, tiny_decay) == 1);
    CHECK(kernel_value("AA", "AA", {KernelType::subsequence, 2, tiny_decay}) == 1);
}

TEST_CASE("Published degree weights and hand-calculated matches are retained",
          "[kernel][characterization]") {
    using namespace stringkernel;
    CHECK(degree_weight(2, 1) == approximately(2.0 / 3));
    CHECK(degree_weight(2, 2) == approximately(1.0 / 3));
    CHECK(weighted_degree("AC", "AC", 2) == approximately(5.0 / 3));
    CHECK(weighted_degree("AC", "AD", 2) == approximately(2.0 / 3));
    CHECK(kernel_value("AC", "AD", {KernelType::weighted_degree, 2}) == approximately(2.0 / 5));
    CHECK(weighted_degree("AAA", "AAB", 2) == approximately(5.0 / 3));
    CHECK(weighted_degree("AAA", "ABA", 2) == approximately(4.0 / 3));
    CHECK(weighted_degree("ABC", "BCA", 3) == 0);
    CHECK(weighted_degree("AC", "ac", 2) == 0);
}

TEST_CASE("Requested degree and unequal-length normalization preserve legacy mathematics",
          "[kernel][characterization]") {
    using namespace stringkernel;
    CHECK(kernel_value("AA", "AB", {KernelType::weighted_degree, 2}) == approximately(2.0 / 5));
    CHECK(kernel_value("AA", "AB", {KernelType::weighted_degree, 5}) == approximately(5.0 / 14));
    CHECK(weighted_degree("A", "AA", 1) == 1);
    CHECK(weighted_degree_self("AA", 1) == 2);
    CHECK(kernel_value("A", "AA", {KernelType::weighted_degree, 1}) ==
          approximately(1 / std::sqrt(2.0)));
}

TEST_CASE("Weighted degree agrees with an independent explicit feature map", "[kernel][oracle]") {
    const auto sequences = short_sequences();
    for (const auto degree : {1, 2, 3, 5}) {
        for (const auto& a : sequences) {
            CAPTURE(degree, a);
            CHECK(stringkernel::weighted_degree_self(a, degree) ==
                  approximately(oracle_raw(a, a, degree)).margin(1e-13));
            for (const auto& b : sequences) {
                CAPTURE(b);
                CHECK(stringkernel::weighted_degree(a, b, degree) ==
                      approximately(oracle_raw(a, b, degree)).margin(1e-13));
                const auto value = stringkernel::kernel_value(
                    a, b, {stringkernel::KernelType::weighted_degree, degree});
                CHECK(value == approximately(oracle_normalized(a, b, degree)).margin(1e-13));
                CHECK(std::isfinite(value));
                CHECK(value >= 0);
                CHECK(value <= 1 + 1e-13);
                CHECK(value ==
                      approximately(stringkernel::kernel_value(
                                        b, a, {stringkernel::KernelType::weighted_degree, degree}))
                          .margin(1e-13));
            }
        }
    }
}

TEST_CASE("Invalid configuration and weight orders are rejected", "[kernel][validation]") {
    using namespace stringkernel;
    for (const auto degree : {0, -1, std::numeric_limits<int>::min()}) {
        CAPTURE(degree);
        CHECK_THROWS_AS((KernelConfig{KernelType::weighted_degree, degree}.validate()),
                        std::invalid_argument);
        CHECK_THROWS_AS(weighted_degree("A", "A", degree), std::invalid_argument);
        CHECK_THROWS_AS(weighted_degree_self("A", degree), std::invalid_argument);
        CHECK_THROWS_AS(degree_weight(degree, 1), std::invalid_argument);
    }
    CHECK_THROWS_AS((KernelConfig{static_cast<KernelType>(99), 5}.validate()),
                    std::invalid_argument);
    CHECK_THROWS_AS(kernel_value("A", "A", {static_cast<KernelType>(99), 5}),
                    std::invalid_argument);
    CHECK_THROWS_AS(degree_weight(5, 0), std::invalid_argument);
    CHECK_THROWS_AS(degree_weight(5, 6), std::invalid_argument);
    CHECK_THROWS_AS(degree_weight(5, std::numeric_limits<std::size_t>::max()),
                    std::invalid_argument);
}

TEST_CASE("Maximum integer degree needs only features present in the inputs",
          "[kernel][large-degree]") {
    using namespace stringkernel;
    const auto degree = std::numeric_limits<int>::max();
    CHECK_NOTHROW((KernelConfig{KernelType::weighted_degree, degree}.validate()));
    CHECK(degree_weight(degree, 1) == approximately(2.0 / (double(degree) + 1)));
    CHECK(degree_weight(degree, static_cast<std::size_t>(degree)) > 0);
    CHECK(weighted_degree("AC", "AC", degree) == approximately(oracle_raw("AC", "AC", degree)));
    CHECK(weighted_degree_self("AC", degree) == approximately(oracle_raw("AC", "AC", degree)));
    CHECK(kernel_value("AC", "AD", {KernelType::weighted_degree, degree}) ==
          approximately(oracle_normalized("AC", "AD", degree)));
}

TEST_CASE("Kernel evaluation has no configuration-dependent or sequence-dependent cache",
          "[kernel][reuse]") {
    using namespace stringkernel;
    KernelConfig config;
    std::string sequence = "AA";
    const auto initial = kernel_value(sequence, "AB", config);
    config.degree = 1;
    CHECK(kernel_value(sequence, "AB", config) == approximately(0.5));
    sequence = "AB";
    CHECK(kernel_value(sequence, "AB", config) == approximately(1));
    config.degree = 5;
    CHECK(kernel_value("AA", "AB", config) == approximately(initial));
}

TEST_CASE("Empty weighted-degree features normalize to zero", "[kernel][normalization]") {
    using namespace stringkernel;
    CHECK(weighted_degree("", "A", 5) == 0);
    CHECK(weighted_degree("A", "", 5) == 0);
    CHECK(weighted_degree_self("", 5) == 0);
    CHECK(kernel_value("", "A", {}) == 0);
    CHECK(kernel_value("", "", {}) == 0);
    CHECK(normalize_kernel(0, 0, 0) == 0);
    CHECK(normalize_kernel(1, 0, 1) == 0);
    CHECK(normalize_kernel(1, 1, 0) == 0);
}

TEST_CASE("Normalization rejects negative and nonfinite inputs", "[kernel][normalization]") {
    using stringkernel::normalize_kernel;
    for (const auto invalid :
         {-1.0, std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::quiet_NaN()}) {
        CHECK_THROWS_AS(normalize_kernel(invalid, 1, 1), std::invalid_argument);
        CHECK_THROWS_AS(normalize_kernel(1, invalid, 1), std::invalid_argument);
        CHECK_THROWS_AS(normalize_kernel(1, 1, invalid), std::invalid_argument);
        CHECK_THROWS_AS(normalize_kernel(invalid, 0, 0), std::invalid_argument);
    }
}

TEST_CASE("Normalization avoids overflow and underflow in self-kernel products",
          "[kernel][normalization]") {
    using stringkernel::normalize_kernel;
    const auto maximum = std::numeric_limits<double>::max();
    const auto minimum = std::numeric_limits<double>::min();
    const auto denormal = std::numeric_limits<double>::denorm_min();
    CHECK(normalize_kernel(maximum, maximum, maximum) == approximately(1));
    CHECK(normalize_kernel(minimum, minimum, minimum) == approximately(1));
    CHECK(normalize_kernel(denormal, denormal, denormal) == approximately(1));
    CHECK(normalize_kernel(1, minimum, maximum) ==
          approximately(1 / (std::sqrt(minimum) * std::sqrt(maximum))));
    CHECK(std::isfinite(normalize_kernel(maximum, minimum, maximum)));
    CHECK(normalize_kernel(0, minimum, maximum) == 0);
}

TEST_CASE("Unequal subnormal self-kernels retain their mathematical ratio",
          "[kernel][normalization]") {
    const auto smallest = std::numeric_limits<double>::denorm_min();
    const auto expected = 1 / std::sqrt(2.0);
    CHECK(stringkernel::normalize_kernel(smallest, smallest, 2 * smallest) ==
          approximately(expected));
    CHECK(stringkernel::normalize_kernel(smallest, 2 * smallest, smallest) ==
          approximately(expected));
}

TEST_CASE("DSK remains equal-length positional similarity", "[kernel][distance]") {
    using namespace stringkernel;
    const KernelConfig config{KernelType::distance, 5};
    CHECK(kernel_value("AC", "AC", config) == 1);
    CHECK(kernel_value("AC", "AD", config) == 0.5);
    CHECK(kernel_value("AC", "CA", config) == 0);
    CHECK(kernel_value("AC", "ac", config) == 0);
    CHECK_THROWS_AS(kernel_value("A", "AA", config), std::invalid_argument);
    CHECK_THROWS_AS(kernel_value("AA", "A", config), std::invalid_argument);
    CHECK_THROWS_AS(kernel_value("", "", config), std::invalid_argument);
}
