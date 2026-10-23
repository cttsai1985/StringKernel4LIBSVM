#include "LibsvmWriter.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace stringkernel {
namespace {
bool is_space(unsigned char character) {
    return character == ' ' || character == '\t' || character == '\n' || character == '\r' ||
           character == '\v' || character == '\f';
}

void validate_samples(const Dataset& data) {
    if (data.empty()) {
        throw std::invalid_argument("kernel datasets must not be empty");
    }
    if (data.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("dataset exceeds LIBSVM's integer index limit");
    }
    for (const auto& sample : data) {
        if (sample.sequence.empty() ||
            std::any_of(sample.sequence.begin(), sample.sequence.end(), is_space)) {
            throw std::invalid_argument("sequences must be nonempty whitespace-free tokens");
        }
        if (!std::isfinite(sample.label)) {
            throw std::invalid_argument("labels must be finite");
        }
    }
}

void validate_lengths(const Dataset& reference, const Dataset& target) {
    const auto length = reference.front().sequence.size();
    for (const auto* data : {&reference, &target}) {
        if (std::any_of(data->begin(), data->end(), [length](const Sample& sample) {
                return sample.sequence.size() != length;
            })) {
            throw std::invalid_argument("DSK requires equal sequence lengths in all datasets");
        }
    }
}

std::vector<double> self_values(const Dataset& data, int degree) {
    std::vector<double> result;
    result.reserve(data.size());
    for (const auto& sample : data) {
        result.push_back(weighted_degree_self(sample.sequence, degree));
    }
    return result;
}

void require_output(const std::ostream& output) {
    if (!output) {
        throw std::runtime_error("failed to write precomputed kernel output");
    }
}
} // namespace

void write_precomputed(std::ostream& output, const Dataset& reference, const Dataset& target,
                       const KernelConfig& config) {
    config.validate();
    validate_samples(reference);
    validate_samples(target);
    if (config.type == KernelType::distance) {
        validate_lengths(reference, target);
    }
    require_output(output);
    const bool weighted = config.type == KernelType::weighted_degree;
    const auto reference_self =
        weighted ? self_values(reference, config.degree) : std::vector<double>{};
    const auto separate_self = weighted && &reference != &target
                                   ? self_values(target, config.degree)
                                   : std::vector<double>{};
    const auto& target_self = &reference == &target ? reference_self : separate_self;
    std::ostringstream row;
    row.imbue(std::locale::classic());
    row << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < target.size(); ++i) {
        row.str("");
        row.clear();
        row << target[i].label << " 0:" << i + 1;
        for (std::size_t j = 0; j < reference.size(); ++j) {
            const double value =
                weighted ? normalize_kernel(weighted_degree(target[i].sequence,
                                                            reference[j].sequence, config.degree),
                                            target_self[i], reference_self[j])
                         : kernel_value(target[i].sequence, reference[j].sequence, config);
            row << ' ' << j + 1 << ':' << value;
        }
        output << row.str() << '\n';
        require_output(output);
    }
    output.flush();
    require_output(output);
}
} // namespace stringkernel
