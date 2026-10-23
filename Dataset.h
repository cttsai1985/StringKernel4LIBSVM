#ifndef STRINGKERNEL_DATASET_H
#define STRINGKERNEL_DATASET_H

#include <filesystem>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

namespace stringkernel {

struct Sample {
    std::string sequence;
    double label;
};

using Dataset = std::vector<Sample>;

Dataset read_dataset(std::istream& input, std::string_view source = "<stream>");
Dataset load_dataset(const std::filesystem::path& path);

} // namespace stringkernel

#endif
