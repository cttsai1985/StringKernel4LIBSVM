#include "Dataset.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

namespace {

class CommaDecimal : public std::numpunct<char> {
  protected:
    char do_decimal_point() const override { return ','; }
};

class BrokenInput : public std::streambuf {
  public:
    explicit BrokenInput(std::string prefix = {}) : prefix_(std::move(prefix)) {
        setg(prefix_.data(), prefix_.data(), prefix_.data() + prefix_.size());
    }

  protected:
    int_type underflow() override { throw std::ios_base::failure("device failed"); }

  private:
    std::string prefix_;
};

} // namespace

TEST_CASE("datasets retain sequence bytes and finite numeric labels", "[dataset]") {
    std::istringstream input("  ACt?\t+1.25\r\n\n \t\r\nXX -2e2\nZ .5");
    const auto result = stringkernel::read_dataset(input);
    REQUIRE(result.size() == 3);
    CHECK(result[0].sequence == "ACt?");
    CHECK(result[0].label == 1.25);
    CHECK(result[1].sequence == "XX");
    CHECK(result[1].label == -200.0);
    CHECK(result[2].sequence == "Z");
    CHECK(result[2].label == 0.5);
}

TEST_CASE("datasets read records longer than the legacy buffer", "[dataset]") {
    const std::string sequence(100000, 'A');
    std::istringstream input(sequence + "\t1\n");
    const auto result = stringkernel::read_dataset(input);
    REQUIRE(result.size() == 1);
    CHECK(result.front().sequence == sequence);
}

TEST_CASE("dataset numbers use the classic locale without modifying the input", "[dataset]") {
    const std::locale comma(std::locale::classic(), new CommaDecimal);
    std::istringstream input("AC 1.25\n");
    input.imbue(comma);
    const auto result = stringkernel::read_dataset(input);
    REQUIRE(result.size() == 1);
    CHECK(result.front().label == 1.25);
    CHECK(std::use_facet<std::numpunct<char>>(input.getloc()).decimal_point() == ',');
}

TEST_CASE("dataset labels retain the finite double range and signed zero", "[dataset]") {
    const std::vector<double> labels = {
        std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::min(), std::numeric_limits<double>::denorm_min(), -0.0};
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto label : labels) {
        text << "A " << label << '\n';
    }
    std::istringstream input(text.str());
    const auto result = stringkernel::read_dataset(input);
    REQUIRE(result.size() == labels.size());
    for (std::size_t index = 0; index < labels.size(); ++index) {
        CHECK(result[index].label == labels[index]);
    }
    CHECK(std::signbit(result.back().label));
}

TEST_CASE("exception-enabled streams allow ordinary EOF", "[dataset]") {
    for (const auto text : {"AC 1\n", "AC 1"}) {
        std::istringstream input(text);
        input.exceptions(std::ios::badbit | std::ios::failbit | std::ios::eofbit);
        const auto result = stringkernel::read_dataset(input);
        REQUIRE(result.size() == 1);
        CHECK(result.front().sequence == "AC");
    }
}

TEST_CASE("malformed dataset records identify the source and physical line", "[dataset]") {
    const std::vector<std::string> records = {
        "AC",     "AC nope", "AC 1 trailing", "AC 1.2suffix", "AC 1,5", "AC nan", "AC NaN",
        "AC inf", "AC -inf", "AC 1e999",      "AC 1e-999",    "AC --1", "AC 0x10"};
    for (const auto& record : records) {
        CAPTURE(record);
        std::istringstream input("A 1\n \t\r\n" + record + "\n");
        REQUIRE_THROWS_WITH(stringkernel::read_dataset(input, "broken data.txt"),
                            Catch::Matchers::ContainsSubstring("broken data.txt:3"));
    }
}

TEST_CASE("empty datasets are rejected", "[dataset]") {
    for (const auto& text : {std::string(), std::string(" \t\n\r\n")}) {
        std::istringstream input(text);
        REQUIRE_THROWS_AS(stringkernel::read_dataset(input), std::runtime_error);
    }
}

TEST_CASE("dataset read failures carry source context", "[dataset]") {
    SECTION("bad streams") {
        std::istringstream input("AC 1\n");
        input.setstate(std::ios::badbit);
        REQUIRE_THROWS_WITH(stringkernel::read_dataset(input, "bad.txt"),
                            Catch::Matchers::ContainsSubstring("bad.txt"));
    }
    SECTION("failed streams without EOF") {
        std::istringstream input("AC 1\n");
        input.setstate(std::ios::failbit);
        REQUIRE_THROWS_WITH(stringkernel::read_dataset(input, "failed.txt"),
                            Catch::Matchers::ContainsSubstring("failed.txt"));
    }
    SECTION("underlying device failure") {
        BrokenInput buffer;
        std::istream input(&buffer);
        REQUIRE_THROWS_WITH(stringkernel::read_dataset(input, "device.txt"),
                            Catch::Matchers::ContainsSubstring("device.txt"));
    }
    SECTION("streams configured to throw") {
        BrokenInput buffer;
        std::istream input(&buffer);
        input.exceptions(std::ios::badbit);
        REQUIRE_THROWS_WITH(stringkernel::read_dataset(input, "exception.txt"),
                            Catch::Matchers::ContainsSubstring("exception.txt"));
    }
    SECTION("device failure after a valid record") {
        BrokenInput buffer("AC 1\n");
        std::istream input(&buffer);
        REQUIRE_THROWS_WITH(stringkernel::read_dataset(input, "partial.txt"),
                            Catch::Matchers::ContainsSubstring("partial.txt:2"));
    }
}

TEST_CASE("file loading names the path in failures", "[dataset]") {
    const auto path = std::filesystem::temp_directory_path() / "stringkernel-nonexistent-parent" /
                      "missing data.txt";
    REQUIRE_THROWS_WITH(stringkernel::load_dataset(path),
                        Catch::Matchers::ContainsSubstring(path.string()));
}

TEST_CASE("file loading reads valid content and diagnoses malformed files", "[dataset]") {
    const auto path = std::filesystem::temp_directory_path() / "stringkernel-dataset-fixture.txt";
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() {
            std::error_code error;
            std::filesystem::remove(path, error);
        }
    } cleanup{path};
    {
        std::ofstream output(path);
        output << "AC 1.25\n";
    }
    const auto result = stringkernel::load_dataset(path);
    REQUIRE(result.size() == 1);
    CHECK(result.front().label == 1.25);
    {
        std::ofstream output(path);
        output << "AC invalid\n";
    }
    REQUIRE_THROWS_WITH(stringkernel::load_dataset(path),
                        Catch::Matchers::ContainsSubstring(path.string() + ":1"));
}
