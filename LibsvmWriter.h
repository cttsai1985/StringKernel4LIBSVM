#ifndef STRINGKERNEL_LIBSVM_WRITER_H
#define STRINGKERNEL_LIBSVM_WRITER_H

#include "Dataset.h"
#include "StringKernels.h"

#include <iosfwd>

namespace stringkernel {

void write_precomputed(std::ostream& output, const Dataset& reference, const Dataset& target,
                       const KernelConfig& config = {});

} // namespace stringkernel

#endif
