#include "CliOptions.h"
#include "Dataset.h"
#include "LibsvmWriter.h"
#include "svm.h"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>

int main(int argc, char** argv) {
    try {
        std::vector<std::string_view> arguments;
        for (int i = 1; i < argc; ++i) {
            arguments.emplace_back(argv[i]);
        }
        const auto options = stringkernel::parse_options(arguments);
        if (options.help) {
            std::cout << stringkernel::usage_text();
        } else if (options.version) {
            std::cout << "StringKernel4LIBSVM " << STRINGKERNEL_VERSION << " (LIBSVM "
                      << LIBSVM_VERSION / 100 << '.' << LIBSVM_VERSION % 100 << ")\n";
        } else {
            const auto training = stringkernel::load_dataset(options.train);
            std::optional<stringkernel::Dataset> target;
            if (options.test) {
                target = stringkernel::load_dataset(*options.test);
            }
            stringkernel::write_precomputed(std::cout, training, target ? *target : training,
                                            options.kernel);
        }
        std::cout.flush();
        if (!std::cout) {
            throw std::runtime_error("failed to write standard output");
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "StringKernelConverter: " << error.what() << "\nUse --help for usage.\n";
        return EXIT_FAILURE;
    }
}
