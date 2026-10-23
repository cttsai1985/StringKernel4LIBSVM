#include "Dataset.h"
#include "LibsvmWriter.h"
#include "StringKernels.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

namespace {
using namespace stringkernel;

std::vector<std::vector<double>> values(const std::string& text) {
    std::vector<std::vector<double>> matrix;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        double label = 0;
        row >> label;
        std::vector<double> entries{label};
        std::string token;
        while (row >> token) {
            entries.push_back(std::stod(token.substr(token.find(':') + 1)));
        }
        matrix.push_back(entries);
    }
    return matrix;
}

class CommaDecimal : public std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
};

class FailingBuffer : public std::streambuf {
    int_type overflow(int_type) override { return traits_type::eof(); }
};

class FailedFlushBuffer : public std::stringbuf {
    int sync() override { return -1; }
};
} // namespace

TEST_CASE("writer pins normalized degree defaults and retains zero columns") {
    const Dataset reference{{"AC", 1}, {"AD", -1}, {"ZZ", 2}};
    std::ostringstream output;
    write_precomputed(output, reference, reference, {KernelType::weighted_degree, 2});
    const auto matrix = values(output.str());
    REQUIRE(matrix.size() == 3);
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        REQUIRE(matrix[i].size() == 5);
        CHECK(matrix[i][1] == static_cast<double>(i + 1));
        CHECK(matrix[i][i + 2] == Catch::Approx(1.0));
    }
    CHECK(matrix[0][3] == Catch::Approx(0.4));
    CHECK(output.str().find(" 3:0\n") != std::string::npos);
    CHECK(matrix[1][2] == Catch::Approx(matrix[0][3]));
}

TEST_CASE("test matrices use reference columns and one-based target serials") {
    const Dataset reference{{"AC", 8}, {"AD", 9}};
    const Dataset target{{"AD", 0.12345678901234566}, {"AC", 4}, {"ZZ", 7}};
    std::ostringstream output;
    write_precomputed(output, reference, target, {KernelType::weighted_degree, 2});
    const auto matrix = values(output.str());
    REQUIRE(matrix.size() == 3);
    CHECK(matrix[0][0] == target[0].label);
    CHECK(matrix[0][2] == Catch::Approx(0.4));
    CHECK(matrix[0][3] == Catch::Approx(1));
    CHECK(matrix[1][2] == Catch::Approx(1));
    CHECK(matrix[2][1] == 3);
    CHECK(matrix[2][2] == 0);
    CHECK(matrix[2][3] == 0);
    CHECK(output.str().find("0:1 1:") != std::string::npos);
    CHECK(output.str().find("0:3 1:0 2:0\n") != std::string::npos);
}

TEST_CASE("permutations move columns and target rows independently") {
    const Dataset reference{{"AA", 1}, {"BB", 2}};
    const Dataset reversed{{"BB", 2}, {"AA", 1}};
    const Dataset target{{"BB", 3}};
    std::ostringstream first, second;
    write_precomputed(first, reference, target);
    write_precomputed(second, reversed, target);
    CHECK(values(first.str())[0][2] == values(second.str())[0][3]);
    CHECK(values(first.str())[0][3] == values(second.str())[0][2]);
}

TEST_CASE("writer recomputes self values for each configuration and dataset") {
    Dataset reference{{"AC", 1}, {"AD", -1}};
    std::ostringstream degree_one, degree_two, changed, single;
    write_precomputed(degree_one, reference, reference, {KernelType::weighted_degree, 1});
    write_precomputed(degree_two, reference, reference, {KernelType::weighted_degree, 2});
    CHECK(values(degree_one.str())[0][3] == Catch::Approx(0.5));
    CHECK(values(degree_two.str())[0][3] == Catch::Approx(0.4));
    reference[1].sequence = "AC";
    write_precomputed(changed, reference, reference);
    CHECK(values(changed.str())[0][3] == Catch::Approx(1));
    write_precomputed(single, reference, Dataset{reference[1]});
    CHECK(values(single.str())[0][2] == values(changed.str())[1][2]);
}

TEST_CASE("writer supports unequal WDSK lengths and equal-length DSK") {
    const Dataset reference{{"AA", 1}, {"A", 2}};
    std::ostringstream wd;
    write_precomputed(wd, reference, reference);
    CHECK(values(wd.str())[0][3] == Catch::Approx(values(wd.str())[1][2]));
    std::ostringstream dsk;
    write_precomputed(dsk, Dataset{{"AB", 1}}, Dataset{{"AC", 2}},
                      {KernelType::distance, 5});
    CHECK(values(dsk.str())[0][2] == Catch::Approx(0.5));
}

TEST_CASE("invalid datasets are rejected before any output") {
    const Dataset good{{"AA", 1}};
    const std::vector<Dataset> invalid{{}, {{"", 1}}, {{"A A", 1}},
                                       {{"AA", std::numeric_limits<double>::infinity()}},
                                       {{"AA", std::numeric_limits<double>::quiet_NaN()}}};
    for (const auto& dataset : invalid) {
        std::ostringstream output;
        CHECK_THROWS_AS(write_precomputed(output, good, dataset), std::invalid_argument);
        CHECK(output.str().empty());
        CHECK_THROWS_AS(write_precomputed(output, dataset, good), std::invalid_argument);
        CHECK(output.str().empty());
    }
    std::ostringstream output;
    CHECK_THROWS_AS(write_precomputed(output, good, Dataset{{"AA", 1}, {"A", 2}},
                                     {KernelType::distance, 5}), std::invalid_argument);
    CHECK(output.str().empty());
    CHECK_THROWS_AS(write_precomputed(output, good, good, {KernelType::weighted_degree, 0}),
                    std::invalid_argument);
    CHECK(output.str().empty());
}

TEST_CASE("writer uses classic locale and round-trip precision without changing the caller") {
    std::ostringstream output;
    output.imbue(std::locale(std::locale::classic(), new CommaDecimal));
    output << std::fixed << std::setprecision(2);
    const auto flags = output.flags();
    const auto precision = output.precision();
    const Dataset data{{"AA", 0.12345678901234566}};
    write_precomputed(output, data, data);
    CHECK(output.str().find(',') == std::string::npos);
    CHECK(values(output.str())[0][0] == data[0].label);
    CHECK(output.flags() == flags);
    CHECK(output.precision() == precision);
    CHECK(std::use_facet<std::numpunct<char>>(output.getloc()).decimal_point() == ',');
}

TEST_CASE("writer reports initially bad streams, write failures, and flush failures") {
    const Dataset data{{"AA", 1}};
    std::ostringstream bad;
    bad.setstate(std::ios::badbit);
    CHECK_THROWS_AS(write_precomputed(bad, data, data), std::runtime_error);
    FailingBuffer buffer;
    std::ostream output(&buffer);
    CHECK_THROWS_AS(write_precomputed(output, data, data), std::runtime_error);
    FailedFlushBuffer flush_buffer;
    std::ostream flushed(&flush_buffer);
    CHECK_THROWS_AS(write_precomputed(flushed, data, data), std::runtime_error);
}
