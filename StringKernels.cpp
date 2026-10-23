#include "StringKernels.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace stringkernel {
namespace {

void validate_degree(int degree) {
    if (degree <= 0) {
        throw std::invalid_argument("kernel degree must be a positive integer");
    }
}

double weight_scale(int degree) {
    const auto d = static_cast<double>(degree);
    return 2.0 / (d * (d + 1.0));
}

double weight_at(int degree, std::size_t order, double scale) {
    return (static_cast<double>(degree) - static_cast<double>(order) + 1.0) * scale;
}

} // namespace

void KernelConfig::validate() const {
    validate_degree(degree);
    if (type != KernelType::distance && type != KernelType::weighted_degree) {
        throw std::invalid_argument("kernel type must be 0 (DSK) or 1 (WDSK)");
    }
}

double degree_weight(int degree, std::size_t order) {
    validate_degree(degree);
    if (order == 0 || order > static_cast<std::size_t>(degree)) {
        throw std::invalid_argument("substring order must be between one and kernel degree");
    }
    return weight_at(degree, order, weight_scale(degree));
}

double weighted_degree(std::string_view a, std::string_view b, int degree) {
    validate_degree(degree);
    const auto length = std::min(a.size(), b.size());
    const auto maximum_order = static_cast<std::size_t>(degree);
    const auto scale = weight_scale(degree);
    double score = 0;

    for (std::size_t position = 0; position < length; ++position) {
        const auto limit = std::min(maximum_order, length - position);
        for (std::size_t order = 1; order <= limit; ++order) {
            if (a[position + order - 1] != b[position + order - 1]) {
                break;
            }
            score += weight_at(degree, order, scale);
        }
    }
    return score;
}

double weighted_degree_self(std::string_view sequence, int degree) {
    validate_degree(degree);
    const auto maximum_order = std::min(sequence.size(), static_cast<std::size_t>(degree));
    const auto scale = weight_scale(degree);
    double score = 0;

    for (std::size_t order = 1; order <= maximum_order; ++order) {
        score += static_cast<double>(sequence.size() - order + 1) * weight_at(degree, order, scale);
    }
    return score;
}

double normalize_kernel(double value, double self_a, double self_b) {
    if (!std::isfinite(value) || !std::isfinite(self_a) || !std::isfinite(self_b) || value < 0 ||
        self_a < 0 || self_b < 0) {
        throw std::invalid_argument("kernel and self-kernel values must be finite and nonnegative");
    }
    if (self_a == 0 || self_b == 0) {
        return 0;
    }
    const auto root_a = std::sqrt(self_a);
    const auto root_b = std::sqrt(self_b);
    const auto denominator = root_a * root_b;
    if (denominator < std::numeric_limits<double>::min()) {
        // Both roots are below one here. Dividing separately retains subnormal
        // ratios without introducing an overflow before a finite result.
        return (value / root_a) / root_b;
    }
    return value / denominator;
}

double kernel_value(std::string_view a, std::string_view b, const KernelConfig& config) {
    config.validate();
    if (config.type == KernelType::weighted_degree) {
        return normalize_kernel(weighted_degree(a, b, config.degree),
                                weighted_degree_self(a, config.degree),
                                weighted_degree_self(b, config.degree));
    }

    if (a.empty() || a.size() != b.size()) {
        throw std::invalid_argument("DSK requires nonempty sequences of equal length");
    }
    std::size_t matches = 0;
    for (std::size_t position = 0; position < a.size(); ++position) {
        if (a[position] == b[position]) {
            ++matches;
        }
    }
    return static_cast<double>(matches) / static_cast<double>(a.size());
}

} // namespace stringkernel
