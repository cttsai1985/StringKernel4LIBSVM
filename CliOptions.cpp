#include "CliOptions.h"

#include <charconv>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>

namespace stringkernel {
namespace {

bool is_value_option(std::string_view option) {
    return option == "-train" || option == "-test" || option == "-k" || option == "-t" ||
           option == "--kernel" || option == "--decay";
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

double parse_decimal(std::string_view token, std::string_view option) {
    double value = 0;
    const auto parsed = std::from_chars(token.data(), token.data() + token.size(), value,
                                        std::chars_format::general);
    if (token.empty() || parsed.ec != std::errc() || parsed.ptr != token.data() + token.size() ||
        !std::isfinite(value)) {
        throw std::invalid_argument(std::string(option) + " requires a finite decimal number");
    }
    return value;
}

KernelType parse_kernel_name(std::string_view value) {
    if (value == "positional") {
        return KernelType::distance;
    }
    if (value == "weighted-degree") {
        return KernelType::weighted_degree;
    }
    if (value == "spectrum") {
        return KernelType::spectrum;
    }
    if (value == "subsequence") {
        return KernelType::subsequence;
    }
    throw std::invalid_argument(
        "--kernel must be positional, weighted-degree, spectrum, or subsequence");
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
        } else if (option == "-t") {
            options.kernel.type = static_cast<KernelType>(parse_integer(value, option));
        } else if (option == "--kernel") {
            options.kernel.type = parse_kernel_name(value);
        } else {
            options.kernel.decay = parse_decimal(value, option);
        }
    }
    if (options.train.empty()) {
        throw std::invalid_argument("-train is required");
    }
    if (seen.count("-t") != 0 && seen.count("--kernel") != 0) {
        throw std::invalid_argument("-t and --kernel select the kernel and cannot be combined");
    }
    if (seen.count("--decay") != 0 && options.kernel.type != KernelType::subsequence) {
        throw std::invalid_argument("--decay is only valid with --kernel subsequence");
    }
    options.kernel.validate();
    return options;
}

std::string_view usage_text() {
    return "Usage: StringKernelConverter -train FILE [-test FILE] [-k K] [-t TYPE | --kernel NAME] "
           "[--decay LAMBDA]\n"
           "       StringKernelConverter --help\n"
           "       StringKernelConverter --version\n"
           "\n"
           "Convert sequence/label records to LIBSVM precomputed kernels.\n"
           "  -train FILE  Training reference sequences and labels (required)\n"
           "  -test FILE   Target sequences and labels (defaults to training data)\n"
           "  -k K         Positive kernel length/order (default: 5)\n"
           "  -t TYPE      Stable numeric kernel ID: 0 positional, 1 weighted degree, 2 spectrum, "
           "3 subsequence\n"
           "  --kernel     positional, weighted-degree, spectrum, or subsequence\n"
           "  --decay      SSK gap decay, finite in (0, 1] (default: 0.5)\n"
           "  --help       Show this help\n"
           "  --version    Show project and bundled LIBSVM versions\n";
}

} // namespace stringkernel
