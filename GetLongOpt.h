#ifndef STRINGKERNEL_OPTIONS_H
#define STRINGKERNEL_OPTIONS_H

#include "StringKernels.h"

#include <filesystem>
#include <optional>
#include <string_view>
#include <vector>

namespace stringkernel {

struct Options {
    std::filesystem::path train;
    std::optional<std::filesystem::path> test;
    KernelConfig kernel;
    bool help = false;
    bool version = false;
};

Options parse_options(const std::vector<std::string_view>& arguments);
std::string_view usage_text();

} // namespace stringkernel

#endif
