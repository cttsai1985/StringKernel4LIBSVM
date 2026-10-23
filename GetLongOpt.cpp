#include "GetLongOpt.h"

#include <charconv>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>

namespace stringkernel {
namespace {

bool is_value_option(std::string_view option) {
    return option == "-train" || option == "-test" || option == "-k" || option == "-t";
}

int parse_integer(std::string_view token, std::string_view option) {
    const auto invalid_integer = [&option]() {
        return std::invalid_argument(std::string(option) +
                                     " requires a complete integer in int range");
    };
    const bool positive_sign = !token.empty() && token.front() == '+';
    if (positive_sign) {
        token.remove_prefix(1);
    }
    if (token.empty() || (positive_sign && token.front() == '-')) {
        throw invalid_integer();
    }
    int value = 0;
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    if (result.ec != std::errc() || result.ptr != token.data() + token.size()) {
        throw invalid_integer();
    }
    return value;
}

std::filesystem::path parse_path(std::string_view value, std::string_view option) {
    if (value.empty() || is_value_option(value) || value == "--help" || value == "--version") {
        throw std::invalid_argument(std::string(option) + " requires a file path");
    }
    return std::filesystem::path(value);
}

} // namespace

Options parse_options(const std::vector<std::string_view>& arguments) {
    Options options;
    if (arguments.size() == 1 && arguments.front() == "--help") {
        options.help = true;
        return options;
    }
    if (arguments.size() == 1 && arguments.front() == "--version") {
        options.version = true;
        return options;
    }

    std::set<std::string_view> seen;
    for (std::size_t index = 0; index < arguments.size(); ++index) {
        const auto option = arguments[index];
        if (option == "--help" || option == "--version") {
            throw std::invalid_argument(std::string(option) + " must be used alone");
        }
        if (!is_value_option(option)) {
            throw std::invalid_argument("unknown option: " + std::string(option));
        }
        if (!seen.insert(option).second) {
            throw std::invalid_argument("duplicate option: " + std::string(option));
        }
        if (++index == arguments.size()) {
            throw std::invalid_argument("missing value for " + std::string(option));
        }
        const auto value = arguments[index];
        if (option == "-train") {
            options.train = parse_path(value, option);
        } else if (option == "-test") {
            options.test = parse_path(value, option);
        } else if (option == "-k") {
            options.kernel.degree = parse_integer(value, option);
        } else {
            options.kernel.type = static_cast<KernelType>(parse_integer(value, option));
        }
    }
    if (options.train.empty()) {
        throw std::invalid_argument("-train is required");
    }
    options.kernel.validate();
    return options;
}

std::string_view usage_text() {
    return "Usage: StringKernelConverter -train FILE [-test FILE] [-k DEGREE] [-t TYPE]\n"
           "       StringKernelConverter --help\n"
           "       StringKernelConverter --version\n"
           "\n"
           "Convert sequence/label records to LIBSVM precomputed kernels.\n"
           "  -train FILE  Training reference sequences and labels (required)\n"
           "  -test FILE   Target sequences and labels (defaults to training data)\n"
           "  -k DEGREE    Positive maximum substring degree (default: 5)\n"
           "  -t TYPE      0: equal-length positional similarity; 1: weighted degree (default)\n"
           "  --help       Show this help\n"
           "  --version    Show project and bundled LIBSVM versions\n";
}

} // namespace stringkernel
