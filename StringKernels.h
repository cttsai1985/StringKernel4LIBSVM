#ifndef STRINGKERNEL_STRING_KERNELS_H
#define STRINGKERNEL_STRING_KERNELS_H

#include <cstddef>
#include <string_view>

namespace stringkernel {

enum class KernelType {
    distance = 0,
    weighted_degree = 1,
    spectrum = 2,
    subsequence = 3,
};

struct KernelConfig {
    KernelType type = KernelType::weighted_degree;
    int degree = 5;
    double decay = 0.5;

    void validate() const;
};

// The requested degree determines weights even when a sequence is shorter.
double degree_weight(int degree, std::size_t order);
double weighted_degree(std::string_view a, std::string_view b, int degree);
double weighted_degree_self(std::string_view sequence, int degree);
double spectrum_kernel(std::string_view a, std::string_view b, int k);
double spectrum_self(std::string_view sequence, int k);
double subsequence_kernel(std::string_view a, std::string_view b, int k, double decay);
double subsequence_self(std::string_view sequence, int k, double decay);

// DSK's raw value is its normalized positional score; its self value is one.
double raw_kernel(std::string_view a, std::string_view b, const KernelConfig& config);
double self_kernel(std::string_view sequence, const KernelConfig& config);

// Nonfinite or negative inputs are invalid; a zero self-kernel yields zero.
double normalize_kernel(double value, double self_a, double self_b);
double kernel_value(std::string_view a, std::string_view b, const KernelConfig& config);

} // namespace stringkernel

#endif
