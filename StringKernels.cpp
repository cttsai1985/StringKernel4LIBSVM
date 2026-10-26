#include "StringKernels.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace stringkernel {
namespace {

void validate_degree(int degree) {
    if (degree <= 0) {
        throw std::invalid_argument("kernel degree must be a positive integer");
    }
}

void validate_decay(double decay) {
    if (!std::isfinite(decay) || decay <= 0 || decay > 1) {
        throw std::invalid_argument("subsequence decay must be finite and in (0, 1]");
    }
}

double weight_scale(int degree) {
    const auto d = static_cast<double>(degree);
    return 2.0 / (d * (d + 1.0));
}

double weight_at(int degree, std::size_t order, double scale) {
    return (static_cast<double>(degree) - static_cast<double>(order) + 1.0) * scale;
}

template <typename T> double checked_result(T value, const char* message) {
    if (!std::isfinite(value) || value < 0 ||
        value > static_cast<T>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<double>(value);
}

using Counts = std::unordered_map<std::string_view, std::size_t>;

Counts spectrum_counts(std::string_view sequence, int k) {
    validate_degree(k);
    Counts result;
    const auto order = static_cast<std::size_t>(k);
    if (order > sequence.size()) {
        return result;
    }
    result.reserve(sequence.size() - order + 1);
    for (std::size_t start = 0; start <= sequence.size() - order; ++start) {
        ++result[sequence.substr(start, order)];
    }
    return result;
}

long double count_dot(const Counts& left, const Counts& right) {
    long double result = 0;
    for (const auto& [token, count] : left) {
        const auto found = right.find(token);
        if (found != right.end()) {
            result += static_cast<long double>(count) * found->second;
        }
    }
    return result;
}

long double subsequence_raw_long(std::string_view a, std::string_view b, int k, double decay) {
    validate_degree(k);
    validate_decay(decay);
    const auto order = static_cast<std::size_t>(k);
    if (order > a.size() || order > b.size()) {
        return 0;
    }
    if (a.size() == std::numeric_limits<std::size_t>::max() ||
        b.size() == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("sequence is too large for subsequence dynamic programming");
    }
    const auto rows = a.size() + 1;
    const auto columns = b.size() + 1;
    if (rows > std::vector<long double>().max_size() / columns) {
        throw std::length_error("subsequence dynamic-programming table is too large");
    }
    const auto cells = rows * columns;
    std::vector<long double> previous(cells, 0);
    std::vector<long double> current(cells, 0);
    for (std::size_t i = 1; i <= a.size(); ++i) {
        for (std::size_t j = 1; j <= b.size(); ++j) {
            if (a[i - 1] == b[j - 1]) {
                previous[i * columns + j] = 1;
            }
        }
    }
    if (order == 1) {
        long double result = 0;
        for (const auto value : previous) {
            result += value;
        }
        return result;
    }

    const auto weight = static_cast<long double>(decay);
    std::vector<long double> gap_row(columns, 0);
    for (std::size_t length = 2; length <= order; ++length) {
        std::fill(current.begin(), current.end(), 0);
        std::fill(gap_row.begin(), gap_row.end(), 0);
        for (std::size_t i = 1; i <= a.size(); ++i) {
            long double prefix = 0;
            for (std::size_t j = 1; j <= b.size(); ++j) {
                prefix = weight * prefix + previous[(i - 1) * columns + (j - 1)];
                gap_row[j] = weight * gap_row[j] + prefix;
                if (a[i - 1] == b[j - 1]) {
                    current[i * columns + j] = gap_row[j];
                }
            }
        }
        previous.swap(current);
    }

    long double result = 0;
    for (const auto value : previous) {
        result += value;
        if (!std::isfinite(result)) {
            throw std::overflow_error("subsequence kernel exceeds numeric range");
        }
    }
    return result;
}

} // namespace

void KernelConfig::validate() const {
    validate_degree(degree);
    switch (type) {
    case KernelType::distance:
    case KernelType::weighted_degree:
    case KernelType::spectrum:
        if (decay != 0.5) {
            throw std::invalid_argument("decay is only configurable for the subsequence kernel");
        }
        break;
    case KernelType::subsequence:
        validate_decay(decay);
        break;
    default:
        throw std::invalid_argument(
            "kernel type must be 0 (DSK), 1 (WDSK), 2 (spectrum), or 3 (subsequence)");
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

double spectrum_kernel(std::string_view a, std::string_view b, int k) {
    return checked_result(count_dot(spectrum_counts(a, k), spectrum_counts(b, k)),
                          "spectrum kernel exceeds numeric range");
}

double spectrum_self(std::string_view sequence, int k) {
    const auto counts = spectrum_counts(sequence, k);
    return checked_result(count_dot(counts, counts), "spectrum self-kernel exceeds numeric range");
}

double subsequence_kernel(std::string_view a, std::string_view b, int k, double decay) {
    return checked_result(subsequence_raw_long(a, b, k, decay),
                          "subsequence kernel exceeds numeric range");
}

double subsequence_self(std::string_view sequence, int k, double decay) {
    return subsequence_kernel(sequence, sequence, k, decay);
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
        return (value / root_a) / root_b;
    }
    return value / denominator;
}

double raw_kernel(std::string_view a, std::string_view b, const KernelConfig& config) {
    config.validate();
    switch (config.type) {
    case KernelType::distance: {
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
    case KernelType::weighted_degree:
        return weighted_degree(a, b, config.degree);
    case KernelType::spectrum:
        return spectrum_kernel(a, b, config.degree);
    case KernelType::subsequence:
        return subsequence_kernel(a, b, config.degree, config.decay);
    }
    throw std::invalid_argument("invalid kernel type");
}

double self_kernel(std::string_view sequence, const KernelConfig& config) {
    config.validate();
    switch (config.type) {
    case KernelType::distance:
        if (sequence.empty()) {
            throw std::invalid_argument("DSK requires nonempty sequences");
        }
        return 1;
    case KernelType::weighted_degree:
        return weighted_degree_self(sequence, config.degree);
    case KernelType::spectrum:
        return spectrum_self(sequence, config.degree);
    case KernelType::subsequence:
        return subsequence_self(sequence, config.degree, config.decay);
    }
    throw std::invalid_argument("invalid kernel type");
}

double kernel_value(std::string_view a, std::string_view b, const KernelConfig& config) {
    config.validate();
    return normalize_kernel(raw_kernel(a, b, config), self_kernel(a, config),
                            self_kernel(b, config));
}

} // namespace stringkernel
