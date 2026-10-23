#ifndef STRINGKERNEL_STRING_KERNELS_H
#define STRINGKERNEL_STRING_KERNELS_H

#include <cstddef>
#include <string_view>

namespace stringkernel {

enum class KernelType { distance = 0, weighted_degree = 1 };

struct KernelConfig {
    KernelType type = KernelType::weighted_degree;
    int degree = 5;

    void validate() const;
};

// The requested degree determines weights even when a sequence is shorter.
double degree_weight(int degree, std::size_t order);
double weighted_degree(std::string_view a, std::string_view b, int degree);
double weighted_degree_self(std::string_view sequence, int degree);

// Nonfinite or negative inputs are invalid; a zero self-kernel yields zero.
double normalize_kernel(double value, double self_a, double self_b);
double kernel_value(std::string_view a, std::string_view b, const KernelConfig& config);

} // namespace stringkernel

#endif
