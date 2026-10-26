#include "CliOptions.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

TEST_CASE("CLI preserves legacy defaults and filenames", "[options]") {
    const auto result = stringkernel::parse_options({"-train", "training data.txt"});
    CHECK(result.train == "training data.txt");
    CHECK_FALSE(result.test.has_value());
    CHECK(result.kernel.type == stringkernel::KernelType::weighted_degree);
    CHECK(result.kernel.degree == 5);
    CHECK(result.kernel.decay == 0.5);
    CHECK_FALSE(result.help);
    CHECK_FALSE(result.version);
}

TEST_CASE("CLI accepts explicit options in any order", "[options]") {
    const auto result = stringkernel::parse_options(
        {"-k", "2", "-test", "test.txt", "-t", "0", "-train", "train.txt"});
    CHECK(result.train == "train.txt");
    REQUIRE(result.test.has_value());
    CHECK(*result.test == "test.txt");
    CHECK(result.kernel.type == stringkernel::KernelType::distance);
    CHECK(result.kernel.degree == 2);
    const auto weighted = stringkernel::parse_options({"-train", "a", "-t", "1"});
    CHECK(weighted.kernel.type == stringkernel::KernelType::weighted_degree);
}

TEST_CASE("CLI supports standalone help and version", "[options]") {
    const auto help = stringkernel::parse_options({"--help"});
    CHECK(help.help);
    CHECK_FALSE(help.version);
    const auto version = stringkernel::parse_options({"--version"});
    CHECK(version.version);
    CHECK_FALSE(version.help);
    CHECK_THAT(std::string(stringkernel::usage_text()),
               Catch::Matchers::ContainsSubstring("-train"));
    CHECK_THAT(std::string(stringkernel::usage_text()),
               Catch::Matchers::ContainsSubstring("-test"));
    CHECK_THAT(std::string(stringkernel::usage_text()), Catch::Matchers::ContainsSubstring("-k"));
    CHECK_THAT(std::string(stringkernel::usage_text()),
               Catch::Matchers::ContainsSubstring("--help"));
    CHECK_THAT(std::string(stringkernel::usage_text()),
               Catch::Matchers::ContainsSubstring("--version"));
}

TEST_CASE("CLI rejects missing unknown duplicate and mixed informational options", "[options]") {
    const std::vector<std::vector<std::string_view>> cases = {
        {},
        {"-test", "a"},
        {"-train"},
        {"-train", ""},
        {"-train", "a", "-test"},
        {"-train", "a", "-test", ""},
        {"-train", "a", "-k"},
        {"-train", "a", "-t"},
        {"-train", "a", "-unknown", "1"},
        {"-train", "a", "extra"},
        {"-train", "a", "-train", "b"},
        {"-train", "a", "-test", "b", "-test", "c"},
        {"-train", "a", "-k", "2", "-k", "3"},
        {"-train", "a", "-t", "0", "-t", "1"},
        {"--help", "--version"},
        {"--help", "--help"},
        {"--help", "-train", "a"},
        {"-train", "a", "--version"},
        {"-train", "-test", "a"}};
    for (const auto& arguments : cases) {
        CAPTURE(arguments);
        REQUIRE_THROWS_AS(stringkernel::parse_options(arguments), std::invalid_argument);
    }
}

TEST_CASE("CLI consumes the complete bounded integer token", "[options]") {
    for (const std::string_view token : {"", "0", "-1", "1.5", "1x", " 2", "2 ", "--1", "++1", "+",
                                         "2147483648", "999999999999999999999999"}) {
        CAPTURE(token);
        REQUIRE_THROWS_AS(stringkernel::parse_options({"-train", "a", "-k", token}),
                          std::invalid_argument);
    }
    for (const std::string_view token : {"-1", "4", "99", "01x", "1.0", "+-0", "2147483648"}) {
        CAPTURE(token);
        REQUIRE_THROWS_AS(stringkernel::parse_options({"-train", "a", "-t", token}),
                          std::invalid_argument);
    }
    const auto largest = std::to_string(std::numeric_limits<int>::max());
    const auto result = stringkernel::parse_options({"-train", "a", "-k", largest});
    CHECK(result.kernel.degree == std::numeric_limits<int>::max());
    REQUIRE_THROWS_AS(stringkernel::parse_options({"-train", "a", "-k", std::string_view{}}),
                      std::invalid_argument);
}

TEST_CASE("CLI parses signed positive integers without numeric prefixes", "[options]") {
    const auto result = stringkernel::parse_options({"-train", "a", "-k", "+2", "-t", "+1"});
    CHECK(result.kernel.degree == 2);
    CHECK(result.kernel.type == stringkernel::KernelType::weighted_degree);
}

TEST_CASE("CLI accepts named and stable numeric kernel selectors", "[options][kernels]") {
    using stringkernel::KernelType;
    CHECK(stringkernel::parse_options({"-train", "a", "--kernel", "positional"}).kernel.type ==
          KernelType::distance);
    CHECK(stringkernel::parse_options({"-train", "a", "--kernel", "weighted-degree"}).kernel.type ==
          KernelType::weighted_degree);
    CHECK(stringkernel::parse_options({"-train", "a", "--kernel", "spectrum"}).kernel.type ==
          KernelType::spectrum);
    const auto subsequence =
        stringkernel::parse_options({"-train", "a", "--kernel", "subsequence"});
    CHECK(subsequence.kernel.type == KernelType::subsequence);
    CHECK(subsequence.kernel.decay == 0.5);
    CHECK(stringkernel::parse_options({"-train", "a", "-t", "0"}).kernel.type ==
          KernelType::distance);
    CHECK(stringkernel::parse_options({"-train", "a", "-t", "1"}).kernel.type ==
          KernelType::weighted_degree);
    CHECK(stringkernel::parse_options({"-train", "a", "-t", "2"}).kernel.type ==
          KernelType::spectrum);
    CHECK(stringkernel::parse_options({"-train", "a", "-t", "3"}).kernel.type ==
          KernelType::subsequence);
    const auto custom = stringkernel::parse_options(
        {"-train", "a", "--kernel", "subsequence", "-k", "3", "--decay", "0.25"});
    CHECK(custom.kernel.degree == 3);
    CHECK(custom.kernel.decay == 0.25);
}

TEST_CASE("CLI restricts decay to subsequence and rejects malformed selectors", "[options]") {
    using stringkernel::parse_options;
    CHECK_THROWS_AS(parse_options({"-train", "a", "--kernel", "spectrum", "-t", "2"}),
                    std::invalid_argument);
    CHECK_THROWS_AS(parse_options({"-train", "a", "--kernel", "spectrum", "--decay", "0.5"}),
                    std::invalid_argument);
    CHECK_THROWS_AS(parse_options({"-train", "a", "--decay", "0.5"}), std::invalid_argument);
    CHECK_THROWS_AS(parse_options({"-train", "a", "--kernel", "unknown"}), std::invalid_argument);
    for (const std::string_view value : {"", "0", "-0.2", "1.1", "nan", "inf", "0.5x"}) {
        CAPTURE(value);
        CHECK_THROWS_AS(parse_options({"-train", "a", "--kernel", "subsequence", "--decay", value}),
                        std::invalid_argument);
    }
    CHECK_THROWS_AS(parse_options({"-train", "a", "--decay", "0.5", "--decay", "0.25"}),
                    std::invalid_argument);
}
