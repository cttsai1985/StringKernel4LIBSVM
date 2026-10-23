# Build, verify, and use the library

## Build and run tests

Use a C++17 compiler, CMake 3.25 or newer, and Python 3.9 or newer. The first test
configuration downloads the pinned Catch2 v3.16.1 archive. After the explicit
LIBSVM bootstrap, production builds with `BUILD_TESTING=OFF` need neither Python
nor downloads.

For Visual Studio 2026, use CMake 4.2 or newer. The CMake 3.25 minimum supports
older toolchains, including Visual Studio 2022.

Run from the repository directory:

```sh
python3 scripts/fetch_libsvm.py
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

In a Windows developer command prompt, use `python scripts/fetch_libsvm.py`
for the bootstrap command. The CMake and CTest commands are the same. Presets use
the generator selected by CMake and provide `Debug` to builds and tests for
generators with multiple configurations. Set `CMAKE_GENERATOR` before the initial
configuration if you need a particular generator.

For a focused regression, list tests with `ctest --preset dev -N`, then run
`ctest --preset dev -R PATTERN`. See [contribution requirements](../CONTRIBUTING.md)
for the expected test-first sequence.

## Bootstrap LIBSVM and work offline

Run `python3 scripts/fetch_libsvm.py` (Windows: `python scripts/fetch_libsvm.py`)
from the repository root. It downloads only the six files in `cmake/libsvm.json`
from the pinned official commit, verifies SHA-256, and stores them in ignored
`third_party/libsvm`. A repeated invocation verifies a complete cache without
network access. An incomplete cache downloads only missing files, publishing
new files only after every required file has passed verification. A bad download
or network failure returns a nonzero status; rerun after resolving connectivity.
For modified cache files, remove the named file and rerun the bootstrap.

For an offline machine, run the bootstrap on a connected machine and transfer the
complete `third_party/libsvm` directory with this checkout. Run the bootstrap again
offline to verify all bytes. CMake also verifies them on every configuration and
prints the setup command for missing or modified files. `--source-dir PATH` lets
you prepare a cache elsewhere; copy that verified directory into
`third_party/libsvm` before configuring. Python is a bootstrap and testing tool;
production configuration, compilation, and installation do not use it.

## Run tests offline

Prepare a verified local copy of the Catch2 v3.16.1 sources. Then configure with
its absolute path:

```sh
cmake --preset dev -DFETCHCONTENT_SOURCE_DIR_CATCH2=/absolute/path/to/Catch2-3.16.1
cmake --build --preset dev
ctest --preset dev
```

The source override bypasses CMake's archive download and digest verification.
Verify the local sources yourself against the pinned archive in
`cmake/Dependencies.cmake`. Prepare the LIBSVM cache as described above.

## Check sanitizers

On a supported Unix platform with GCC or Clang:

```sh
cmake --preset sanitizers
cmake --build --preset sanitizers
ctest --preset sanitizers
```

This configuration instruments project code and its C++ tests with AddressSanitizer
and UndefinedBehaviorSanitizer. It does not alter downloaded LIBSVM. Windows and
unsupported compilers reject the sanitizer option.

## Check coverage

On Linux with GCC, install the same `gcovr` version used in CI and run:

```sh
python3 -m pip install gcovr==8.3
cmake --preset coverage
cmake --build --preset coverage
ctest --preset coverage
python3 -c "from pathlib import Path; Path('coverage').mkdir(exist_ok=True)"
gcovr --root . --filter '(StringKernels|Dataset|LibsvmWriter|CliOptions|StringKernelConverter)\.cpp$' --print-summary --fail-under-line 98 --html-details coverage/index.html --xml coverage/coverage.xml
```

The filter includes project implementation and excludes tests and upstream code.
The line-coverage floor is a regression gate. Consult
[validation evidence](VALIDATION.md) for the measured baseline. Coverage requires
GCC, and the build rejects combined coverage and sanitizer instrumentation.

## Format project C++

Use the pinned formatter version from CI:

```sh
python3 -m pip install clang-format==19.1.7
clang-format --dry-run --Werror Dataset.h Dataset.cpp CliOptions.h CliOptions.cpp StringKernels.h StringKernels.cpp LibsvmWriter.h LibsvmWriter.cpp StringKernelConverter.cpp tests/test_dataset.cpp tests/test_options.cpp tests/test_kernels.cpp tests/test_writer.cpp
```

To apply formatting, replace `--dry-run --Werror` with `-i`. Keep downloaded LIBSVM
and the historical analyzer unchanged.

## Build and install without test dependencies

```sh
python3 scripts/fetch_libsvm.py
cmake -S . -B build/release -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --config Release
cmake --install build/release --config Release --prefix "$PWD/build/install"
```

On Windows, use `--prefix "C:/path/to/install"` instead. Installation follows
CMake's `GNUInstallDirs` variables:

| Location | Contents |
| --- | --- |
| `${CMAKE_INSTALL_BINDIR}` | `StringKernelConverter`, `svm-train`, and `svm-predict`. |
| `${CMAKE_INSTALL_LIBDIR}` | Static `stringkernel` library. |
| `${CMAKE_INSTALL_INCLUDEDIR}/StringKernel4LIBSVM` | The supported public headers. |
| `${CMAKE_INSTALL_LIBDIR}/cmake/StringKernel4LIBSVM` | The `stringkernelTargets.cmake` export and configuration-specific target files. |
| `${CMAKE_INSTALL_DATADIR}/StringKernel4LIBSVM` | `THIRD_PARTY_NOTICES.md` and LIBSVM `COPYRIGHT`. |

The usual directory defaults are `bin`, `lib`, `include`, and `share`. Some
platforms choose a different library directory. No project package configuration
file is installed, so `find_package(StringKernel4LIBSVM)` is not provided.

## Use the C++ library

Link the `stringkernel` CMake target within this project. For an installed tree,
include its `stringkernelTargets.cmake` and link `stringkernel::stringkernel`.
For example, after the default installation above:

```cmake
include("/absolute/path/to/build/install/lib/cmake/StringKernel4LIBSVM/stringkernelTargets.cmake")
target_link_libraries(your_target PRIVATE stringkernel::stringkernel)
```

Adjust `lib` if your install uses another `CMAKE_INSTALL_LIBDIR`.

```cpp
#include "Dataset.h"
#include "LibsvmWriter.h"
#include "StringKernels.h"

#include <iostream>

int main() {
	const stringkernel::Dataset training{{"AC", 1.0}, {"AD", -1.0}};
	const stringkernel::Dataset target{{"AC", 0.0}};
	const stringkernel::KernelConfig config{stringkernel::KernelType::weighted_degree, 2};
	stringkernel::write_precomputed(std::cout, training, target, config);
}
```

`weighted_degree` returns an unnormalized pair value.
`weighted_degree_self` returns its self-kernel. `kernel_value` returns the
normalized configured kernel, and `write_precomputed` streams a complete matrix.
`read_dataset` parses a stream, while `load_dataset` opens and parses a path.
All dataset ownership stays in values and standard containers. The library
neither reads the clock nor changes the process locale.

The writer validates all records before its first row, precomputes self-kernels,
and retains only the current output row. It does not retain the complete matrix.
When the same `Dataset` object is both reference and target, its self-kernels are
reused within that call. Later calls with different configurations calculate new
self-kernels.

## Inspect the verification gates

`ctest` runs independent numerical checks, dataset and option tests, matrix writer
tests, process-level CLI checks, dependency pins, and real LIBSVM workflows.
The integration checks exercise the supplied regression example, a separate
classification fixture, model save and reload, prediction order, and installation.

CI runs the development workflow on Linux, macOS, and Windows. Linux jobs add
sanitizers, formatting, and project coverage. The production build explicitly
disables tests to check the dependency boundary. See
[the workflow](../.github/workflows/ci.yml) for the exact CI commands.
