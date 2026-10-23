# StringKernel4LIBSVM

StringKernel4LIBSVM converts labeled sequences into normalized string-kernel
matrices for LIBSVM. It includes the converter, a C++17 kernel library, and the
unmodified LIBSVM 3.37 training and prediction programs.

The default is the weighted-degree string kernel with degree `5`. Matches depend
on position, so align your sequences before conversion when your application
requires alignment. See [kernel definitions](docs/KERNELS.md).

## Build and run the supplied regression example

Use a C++17 compiler and CMake 3.25 or newer. This production build uses vendored
LIBSVM sources and needs neither Python nor dependency downloads.

On Linux or macOS, run these commands from the repository directory:

```sh
cmake -S . -B build/release -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --config Release
mkdir -p build/example
./build/release/StringKernelConverter -train data.txt -k 5 > build/example/train.kernel
./build/release/StringKernelConverter -train data.txt -test test_data.txt -k 5 > build/example/test.kernel
./build/release/svm-train -s 3 -t 4 build/example/train.kernel build/example/model
./build/release/svm-predict build/example/test.kernel build/example/model build/example/predictions.txt
```

The last command writes one prediction per test record to
`build/example/predictions.txt`. The supplied labels are continuous, so the example
uses epsilon-SVR with `-s 3`. LIBSVM's `-t 4` selects a precomputed kernel.

On Windows, open a Visual Studio developer command prompt with a C++17 toolchain.
For Visual Studio 2026, use CMake 4.2 or newer.

```bat
cmake -S . -B build/release -DBUILD_TESTING=OFF
cmake --build build/release --config Release
sample.bat build\release\Release
```

`sample.bat` writes the same artifacts under `build\example` and stops on errors.
For a generator that puts executables directly in `build\release`, pass that
directory instead. In PowerShell, invoke the batch file as `./sample.bat`. Use
this batch file or `cmd /c` for native output redirection on older PowerShell
versions. Their default `>` encoding can produce UTF-16 text that LIBSVM cannot
read.

## Convert your sequences

Each nonblank line contains a sequence token and a finite numeric label:

```text
ACDE 1
ACDF -1
```

Separate fields with spaces or tabs. Sequences are case-sensitive byte strings.
The parser accepts LF, CRLF, and a final line without a newline. It rejects empty
datasets, missing labels, extra fields, and invalid numbers with a source and
line diagnostic. Save text without a byte-order mark.

Use the same training file, record order, kernel type, and degree for both
conversions. A test matrix always has one column per original training record.
Keep all output columns, including zero values and the `0:serial` field. Do not
run `svm-scale` on these matrices as though they were ordinary feature vectors.

`StringKernelConverter --help` lists the options. The converter's `-t 0` selects
equal-length positional similarity and `-t 1` selects weighted degree. These type
numbers belong to the converter. The bundled `svm-train` needs `-t 4` for either
converted matrix.

## Develop and install

For tests, TDD, sanitizers, coverage, offline test dependencies, installation, and
the library API, see [development instructions](docs/DEVELOPMENT.md) and
[contribution requirements](CONTRIBUTING.md).

The modernized C++ API and stricter input rules are intentional changes. See
[compatibility](docs/COMPATIBILITY.md), [design](docs/DESIGN.md), and
[validation evidence](docs/VALIDATION.md). The incomplete `svmsk_analyser.cpp`
remains unsupported historical source and is excluded from the build.

LIBSVM's license and exact source pin are in
[third-party notices](THIRD_PARTY_NOTICES.md). This modernization does not assign
a license to the historical project.
