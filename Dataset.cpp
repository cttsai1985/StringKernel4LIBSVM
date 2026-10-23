#include "Dataset.h"

#include <charconv>
#include <cmath>
#include <fstream>
#include <istream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace stringkernel {
namespace {

std::runtime_error input_error(std::string_view source, std::size_t line, std::string_view reason) {
    return std::runtime_error(std::string(source) + ":" + std::to_string(line) + ": " +
                              std::string(reason));
}

bool read_line(std::istream& input, std::string& line, std::string_view source,
               std::size_t next_line) {
    line.clear();
    try {
        if (std::getline(input, line)) {
            return true;
        }
    } catch (const std::exception&) {
        if (input.eof() && !input.bad()) {
            return !line.empty();
        }
        throw input_error(source, next_line, "failed to read dataset");
    }
    if (input.bad() || !input.eof()) {
        throw input_error(source, next_line, "failed to read dataset");
    }
    return false;
}

double parse_label(const std::string& token, std::string_view source, std::size_t line) {
    std::string_view number(token);
    if (!number.empty() && number.front() == '+') {
        number.remove_prefix(1);
        if (number.empty() || number.front() == '+' || number.front() == '-') {
            throw input_error(source, line, "label must be a finite decimal double");
        }
    }
    double label = 0;
    const auto result = std::from_chars(number.data(), number.data() + number.size(), label,
                                        std::chars_format::general);
    if (result.ec != std::errc() || result.ptr != number.data() + number.size() ||
        !std::isfinite(label)) {
        throw input_error(source, line, "label must be a finite decimal double");
    }
    return label;
}

} // namespace

Dataset read_dataset(std::istream& input, std::string_view source) {
    Dataset samples;
    std::string line;
    std::size_t line_number = 0;
    while (read_line(input, line, source, line_number + 1)) {
        ++line_number;
        std::istringstream record(line);
        record.imbue(std::locale::classic());
        std::string sequence;
        if (!(record >> sequence)) {
            continue;
        }
        std::string label;
        std::string extra;
        if (!(record >> label) || (record >> extra)) {
            throw input_error(source, line_number, "expected exactly a sequence and label");
        }
        samples.push_back({std::move(sequence), parse_label(label, source, line_number)});
    }
    if (samples.empty()) {
        throw std::runtime_error(std::string(source) + ": dataset is empty");
    }
    return samples;
}

Dataset load_dataset(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error(path.string() + ": could not open dataset");
    }
    return read_dataset(input, path.string());
}

} // namespace stringkernel
