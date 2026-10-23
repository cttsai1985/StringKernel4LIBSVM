# Validation evidence

Tracking issue: [#1](https://github.com/cttsai1985/StringKernel4LIBSVM/issues/1).
This record distinguishes observed results from planned checks.

## RED: process regressions before implementation

Built the original `origin/master` sources in a scratch directory with:

```sh
g++ -std=c++17 -w StringKernelConverter.cpp StringKernels.cpp GetLongOpt.cpp -o StringKernelConverter
python3 tests/test_cli.py --converter /path/to/original/StringKernelConverter
```

Observed failures include successful exit on empty/missing input, acceptance of
invalid degrees/types/options, ignored malformed/nonfinite records, regression
label precision loss, asymmetric unequal-length DSK, and successful exit when
redirected output fails. These tests were written and run before replacing the
converter or implementing the writer. Original valid kernel/CLI characterization
remains a separate compatibility check using parsed values.

Writer API RED command:

```sh
g++ -std=c++17 -I. -fsyntax-only tests/test_writer.cpp
```

Observed compiler error: `LibsvmWriter.h: No such file or directory`. Tests existed
before the new writer API and implementation. Semantic process regressions above
provide failure evidence beyond missing declarations.

Further unit, integration, sanitizer, coverage, and hosted CI results will be
recorded after they actually run.

## Numerical characterization and GREEN results

MEASURED on the original source through a temporary API adapter: the existing
weighted-degree behavior passed the independent feature-map/characterization
suite with 5,474 assertions across three cases. Legacy normalization and DSK
regressions failed 22 assertions. The corrected kernel suite subsequently passed
5,544 assertions across 12 cases under ASan/UBSan and warnings as errors. A
separate unequal-subnormal normalization regression failed both assertions
before its targeted fix and passed afterward.

MEASURED parser suite: 109 assertions across 16 cases passed. Besides initial
missing-API REDs, focused semantic REDs demonstrated silent underflow of a
nonzero label and erroneous acceptance of the malformed integer `+-0`, followed
by GREEN fixes.

MEASURED CLI baseline: 14 test methods produced 29 failing assertions/subtests
against the original executable. Those process contracts subsequently passed
against the modern executable.

An actual old/new conversion comparison used both supplied data files with
default WDSK degree. Parsed values, labels, serials, and matrix entries agreed
within the legacy decimal-rounding bound. MEASURED training output: eight rows,
80 compared values, maximum absolute difference `4.900000014629313e-7`.
MEASURED target output: six rows, 60 compared values, maximum absolute difference
`4.78260869574143e-7`. The comparison asserted a maximum difference no greater
than `5.0000001e-7`; it did not require identical decimal text.

Dependency tests also rejected an isolated copy in which both `svm.cpp` and its
manifest checksum were changed, demonstrating that source/manifest drift cannot
bypass the independently pinned expected hashes.

## Full local verification

MEASURED development, AddressSanitizer/UndefinedBehaviorSanitizer, and GCC coverage
configurations each passed all 39 CTest checks. The suite includes the real
classification and supplied regression workflows, independent matrix/prediction
oracles, model persistence/reload, installation, dependency pins, and CLI errors.
Project sources build with warnings treated as errors.

MEASURED project-only coverage: 260/262 lines (99.2%), 28/28 functions (100%), and
382/514 branches (74.3%). CI enforces a line-coverage floor of 98%; vendored solver
and test code are excluded. Reports are generated under `coverage/` and uploaded
by hosted CI. The unreachable-in-small-fixtures record-count guard remains outside
the exercised input range.

MEASURED documentation checks: the literal README Unix release workflow generated
eight training rows, six target rows, and six finite predictions. A separate
installed-library consumer compiled the documented C++ example through the
exported CMake target and produced the analytical degree-two kernel value.
These are correctness checks, not biological model-quality claims.

Production configuration/build also succeeded with testing disabled, Python
discovery explicitly disabled, and FetchContent disconnected. This configuration
creates no test dependency directory.

In this workspace, Catch2 discovery initially targeted an unwritable inherited
`XDG_RUNTIME_DIR`. The diagnostic was traced to Catch2's temp-file implementation;
local checks use `env -u XDG_RUNTIME_DIR TMPDIR=/tmp ctest --preset dev` (and the
corresponding sanitizer/coverage preset). No product workaround was added.

## Hosted CI and release gates

The first hosted Windows configure exposed a real toolchain mismatch: CMake
3.31.6 predates Visual Studio 2026 support and selected unavailable `nmake` on the
current runner. CI now pins CMake 4.4.4, which supports that generator. Vendor
newline conversion is disabled so independent byte hashes remain valid under
Windows Git settings. Windows CI additionally executes `sample.bat` itself.

Current-head hosted checks and an independent reviewer verdict are required
before merge. Consult the linked PR for their exact SHA, logs, and final result.

## Cross-platform numeric parsing regression

[Hosted run 39899594757](https://github.com/cttsai1985/StringKernel4LIBSVM/actions/runs/39899594757)
established additional RED evidence after the build corrections. macOS rejected
the valid subnormal in `dataset labels retain the finite double range and signed
zero`; macOS and Windows both accepted the hexadecimal token in `malformed
dataset records identify the source and physical line`. Linux passed those
same existing assertions.

The correction uses decimal `std::from_chars` with complete-token, range, and
finite-value checks, preserving a single optional leading plus. It removes
standard-stream parsing differences without relaxing the expected values.
Additional process checks cover subnormal/maximal finite labels and malformed
hexadecimal/double-sign inputs. MEASURED corrected local development suite:
39/39 CTest checks passed. The hosted matrix remains the cross-platform gate.


## Dependency bootstrap and layout cleanup (#3)

Tracking issue: [#3](https://github.com/cttsai1985/StringKernel4LIBSVM/issues/3).
The modernization evidence above records historical paths and remains unchanged.
Current examples are `examples/regression/train_sequences.tsv` and
`examples/regression/test_sequences.tsv`; historical source is
`legacy/paired_sequence_analyzer.cpp`. The active internal option parser is
`CliOptions.h/.cpp`. Independent tests pin the original bytes of both examples
and the historical analyzer. No kernel or solver mathematics changed.

RED: `python3 tests/test_fetch_libsvm.py -v` failed before implementation because
`scripts/fetch_libsvm.py` did not exist. An additional external-response test
raised `http.client.IncompleteRead` through the initial implementation; it passed
after interrupted HTTP responses became actionable bootstrap errors.

MEASURED bootstrap tests: all 19 standard-library tests passed, covering real
byte digests, temporary files, exact upstream pins, complete offline reuse,
partial cache recovery, modified cache rejection, bad hashes, interrupted/network
failures, manifest validation, and CLI diagnostics. HTTP is mocked only at its
external boundary; unit tests make no network requests.

MEASURED CMake preflight tests accept the verified cache and reject missing or
modified sources with the setup command. Production configure/build succeeded
with `BUILD_TESTING=OFF`, Python discovery disabled, and FetchContent disconnected.
Production configuration invokes no dependency fetcher. The formatting check
passed using clang-format 19.1.7.

MEASURED full local verification: development, ASan/UBSan, and GCC coverage each
passed all 40 CTest checks, including the 19-case bootstrap driver, independent
matrix and prediction oracles, actual trainer/predictor workflows, installation,
and dependency/layout pins. Project-only coverage was 260/262 lines (99.2%),
28/28 functions (100%), and 382/508 branches (75.2%); the 98% line gate passed.
The local Catch2 source override and temporary-directory environment remain as
documented above.
