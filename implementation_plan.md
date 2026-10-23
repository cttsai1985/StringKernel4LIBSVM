# StringKernel4LIBSVM modernization plan

Status: authorized by the owner to deepen the design, open a tracking issue,
implement, verify, merge, and close out. See `docs/DESIGN.md` for detailed contracts.

## Goal and working assumptions

Modernize the working sequence-to-kernel converter, integrate the latest stable
LIBSVM, and make numerical behavior, error handling, builds, and documentation
maintainable through test-driven development.

Keep C++ as the implementation language. Interpret “with AI” as AI-assisted
design, implementation, and review. Proposed product scope is the converter,
reusable kernel library, and official LIBSVM training/prediction executables.
Python bindings, new AI product features, and reconstruction of the historical
analyzer would require separate designs.

Risk assessment: high. Numerical output, input validation, ownership, builds,
and dependency distribution all change. Use the Design Solution and TDD skills;
The owner's instruction to refactor and close out authorizes implementation of
this plan; no further planning approval is pending.

## Findings and evidence

These findings come from source inspection and upstream metadata. No builds,
tests, benchmarks, or prediction experiments have been run during planning.

- The repository has no alternate remote branches, open pull requests, build
  system, test suite, or CI configuration. Its default branch is `master`.
- There is no bundled LIBSVM to upgrade. The existing converter already emits
  LIBSVM precomputed-kernel rows. Integration should use this supported boundary.
- Upstream currently identifies **LIBSVM 3.37**, tag `v337`, commit
  `6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb`, as its stable release.
- The weighted-degree formula is already correct for the published equal-length
  definition. The implementation also accepts unequal lengths by considering
  shared positions and normalizing with each full self-kernel.
- `LoadDataset` mismatches array allocation/deallocation, leaves files open,
  uses a fixed-size record buffer, and silently skips malformed records.
  CLI parsing can read beyond `argv`; missing files and unsupported kernel
  selections can report success. Ownership and cached self-kernels are unsafe.
- DSK divides positional matches by the reference sequence length. Unequal
  lengths can therefore produce an asymmetric training matrix.
- `svmsk_analyser.cpp` is incomplete: `OrthogonalArray.h` is absent, paired
  sequence fields and constructor signatures disagree with the supplied
  header, and it requires customized LIBSVM application internals.
- The project declares no license. Retain LIBSVM's copyright and license;
  record project license selection as an owner decision.

Primary references:

- [Official LIBSVM release page](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
- [Pinned LIBSVM source](https://github.com/cjlin1/libsvm/tree/v337)
- [Precomputed-kernel contract](https://github.com/cjlin1/libsvm/blob/v337/README#precomputed-kernels)
- [Public LIBSVM API](https://github.com/cjlin1/libsvm/blob/v337/svm.h)
- [LIBSVM copyright](https://github.com/cjlin1/libsvm/blob/v337/COPYRIGHT)
- [Weighted-degree reference, section 4.2.3](https://www.mlsec.org/docs/2007-lskm.pdf)
- [Catch2 release v3.16.1](https://github.com/catchorg/Catch2/releases/tag/v3.16.1)

## Proposed architecture and contracts

Use C++17 and CMake with CTest. Keep the existing source filenames to limit
unnecessary movement. Build a project kernel library, `StringKernelConverter`,
and the unmodified upstream `svm-train` and `svm-predict` programs.

Responsibilities:

- `StringKernels`: pure pairwise kernel evaluation and validated configuration.
- `Dataset`: value-owned sequences/labels and strict, line-aware input parsing.
- `LibsvmWriter`: streaming matrix output with fixed training-column order.
- `GetLongOpt`: validated CLI options with useful diagnostics.
- `StringKernelConverter`: orchestration and exit-status handling.
- `third_party/libsvm`: unmodified upstream solver sources plus provenance.

Data flow: parse options -> load and validate training/target records -> prepare
self-kernel values -> emit rows -> official `svm-train -t 4` -> `svm-predict`.
Exceptions/errors become stderr diagnostics and nonzero CLI exits. Successful
stdout contains only machine-readable matrix data. Validate input before output
begins, and report output-stream failures.

Compatibility decisions proposed for approval:

- Preserve executable naming, `-train`, `-test`, `-k`, and `-t`; retain default
  WDSK and degree `5`. Add `--help` and `--version`.
- Preserve whitespace-delimited `sequence label` records, case-sensitive byte
  matching, and support for arbitrary nonwhitespace sequence tokens. Document
  LF/CRLF handling; ignore whitespace-only lines. Reject malformed nonempty
  records with a line number, empty datasets, missing values, and nonfinite labels.
- Preserve WDSK normalization and its unequal-length extension. Keep weights
  based on the requested degree even when sequences are shorter; validate a
  positive representable degree and avoid integer-overflow intermediates. Compute
  weights lazily only for substring lengths present in the input; never allocate
  or loop in proportion to an arbitrarily large requested degree.
- Preserve DSK on equal-length sequences. Reject unequal lengths in DSK across
  both training and target datasets; document this intentional validation change.
- Emit all training columns, including zeros, with one-based indices and valid
  `0:<serial>` entries. Test matrices always use the original training order.
- Use locale-independent, round-trip double precision. Parsed values must agree
  with valid legacy output within its rounding tolerance; byte-identical output
  is not promised because legacy output rounded to six decimal places.
- Replace the unsafe internal C++ API with value-owned types and scoped enums.
  Document the internal API change; CLI/data compatibility is the migration target.
- Preserve `svmsk_analyser.cpp` unchanged as unsupported historical source and
  exclude it explicitly from build targets. Document what restoration requires.
  The dormant RBF-over-similarity transformation is outside the supported CLI.

Use vectors, strings, const references, and RAII; remove manual allocation and
configuration-dependent mutable caches. Stream matrix rows instead of retaining
the complete matrix. Preserve numerical behavior before considering optimization;
any future performance claim must be labeled MEASURED or ESTIMATED with evidence.

Vendor the minimal LIBSVM source set, upstream README, and COPYRIGHT at the exact
commit above. Record source URLs and per-file SHA-256 digests and verify them in
CI. This keeps production builds independent of a runtime package manager or
network fetch. Pin Catch2 v3.16.1 and its archive digest for test-only CMake
FetchContent, with a documented local-source override for offline testing.
Python's standard library may drive CLI/integration/provenance tests; Python is
a development dependency only. Discover Python and configure Catch2 only when
`BUILD_TESTING` is enabled. Explicit offline production build:

```sh
cmake -S . -B build -DBUILD_TESTING=OFF
cmake --build build --config Release
```

## Implementation sequence after approval

### Establish build and independent behavior checks

Add explicit CMake targets and a CTest harness before refactoring production
logic. Capture valid legacy sample behavior as characterization evidence; add
hand-calculated cases and an independent test-only feature-map oracle keyed by
position, substring length, and substring, with weighted dot products. Keep
characterization tests distinct from
regressions expected to fail against the old code.

For each defect or changed contract: write its test, observe the expected RED
failure, implement the smallest fix, observe GREEN, then refactor with tests
passing. Commit the sequence so failures and their fixes are reviewable.

ESTIMATED complexity: medium. Main risk: accidentally treating legacy bugs or
rounded output as the mathematical specification.

### Integrate the pinned upstream solver

Add the unmodified vendored sources and project-owned CMake glue. Verify version,
provenance, and digests. Build upstream `.c` CLI sources as required by its C++
build conventions. Retain upstream notices in source and installed artifacts.

ESTIMATED complexity: low to medium. Depends on the build harness; main risk is
platform differences in compiler and installation behavior.

### Refactor under regression tests

Separate configuration, data, mathematics, formatting, and orchestration. Replace
unsafe ownership/parsing, enforce the approved input/length policies, and improve
precision. Retain the supported legacy command line and sample files. Keep
project warnings separate from vendored-source warnings.

ESTIMATED complexity: medium to high. Depends on characterization coverage;
main risks are altered normalization, training-column order, and error behavior.

### Prove training and prediction

Test conversion -> actual upstream training -> persisted model reload -> actual
prediction. Use the supplied continuous-label data for an epsilon-SVR example
(`-s 3 -t 4`) and a separate deterministic classification fixture with expected
held-out labels. Assert finite predictions, dimensions, ordering, and expected
behavior; use explicit numerical tolerances across platforms.

ESTIMATED complexity: medium. Depends on the refactor and solver integration.
These are integration checks, not evidence of biological generalization quality.

### Document, review, and land

Write executable build/use examples for Linux, macOS, and Windows, including
PowerShell encoding guidance. Explain the equation, positional alignment,
normalization, degree and length policies, matrix dimensions, train/test ordering,
and why precomputed matrices must not be scaled as ordinary feature vectors.
Document TDD, dependency updates, numerical tolerances, compatibility changes,
and analyzer limitations.

Review the diff, run the required CI gates, and merge the implementation PR to
`master` only after its current commit passes. Remove its worktree in the same
session after merge. Use incremental commits within the PR rather than leaving
parallel unfinished branches. No data/artifact copies or per-task caches.

ESTIMATED complexity: medium. Depends on the supported commands and behavior
being stable. No elapsed-time or speedup estimate has been measured.

## Verification gates to implement

- Numerical: independent oracle and hand calculations; weight/default pins;
  identical/disjoint/duplicate sequences; normalization, symmetry and diagonal;
  degree boundaries and overflow; unequal lengths under each approved policy.
- Data/CLI: missing files/arguments; unknown options; malformed and nonfinite
  labels; empty input; LF/CRLF; long records; stderr/stdout separation; nonzero
  failures; output precision and write failures; default and version pins.
- Matrix: train/train and test/train dimensions; full ordered columns and zero
  entries; serials; reference-order permutations; batch/single-row consistency.
- LIBSVM: exact source integrity and version pins; real classification and
  regression training, model persistence/reload, and held-out predictions.
- CI: Linux/GCC, macOS/Clang, Windows/MSVC; one Linux AddressSanitizer and
  UndefinedBehaviorSanitizer configuration; formatting and project warnings.
  Publish coverage for project code and inspect uncovered critical/error paths.
  Enforce a regression coverage floor established from the resulting baseline;
  do not count vendored solver code toward project coverage.
- Documentation: run the documented build/convert/train/predict workflow in CI;
  verify installed executables and required third-party notices.

Proposed common local commands after implementation:

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev --output-on-failure
cmake --preset sanitizers
cmake --build --preset sanitizers
ctest --preset sanitizers --output-on-failure
```

## Exact proposed file scope

Modify existing files:

```text
.gitignore
GetLongOpt.h
GetLongOpt.cpp
StringKernels.h
StringKernels.cpp
StringKernelConverter.cpp
README.md
sample.bat
```

Create project/build/documentation files:

```text
implementation_plan.md
AGENTS.md
.clang-format
CMakeLists.txt
CMakePresets.json
cmake/Dependencies.cmake
Dataset.h
Dataset.cpp
LibsvmWriter.h
LibsvmWriter.cpp
.github/workflows/ci.yml
CONTRIBUTING.md
THIRD_PARTY_NOTICES.md
docs/COMPATIBILITY.md
docs/KERNELS.md
docs/DEVELOPMENT.md
docs/DESIGN.md
docs/VALIDATION.md
```

Create tests and dependency files:

```text
tests/CMakeLists.txt
tests/test_kernels.cpp
tests/test_dataset.cpp
tests/test_options.cpp
tests/test_writer.cpp
tests/test_cli.py
tests/test_libsvm_integration.py
tests/test_dependency_pins.py
tests/fixtures/classification_train.txt
tests/fixtures/classification_test.txt
third_party/libsvm/svm.cpp
third_party/libsvm/svm.h
third_party/libsvm/svm-train.c
third_party/libsvm/svm-predict.c
third_party/libsvm/README
third_party/libsvm/COPYRIGHT
third_party/libsvm/UPSTREAM.json
```

Preserve `data.txt`, `test_data.txt`, and `svmsk_analyser.cpp` unchanged. If new
evidence requires files outside this list, update the plan and notify the owner
before editing them. Project license selection is an explicit owner decision;
adding a new project LICENSE is outside this proposed scope.

## Definition of done

- Approved converter contracts are implemented and documented.
- Current pinned upstream LIBSVM builds and performs the documented workflow.
- RED/GREEN evidence exists for meaningful fixes; independent numerical tests
  and real train/predict checks pass on the supported platforms.
- CI gates fail on regressions, dependency drift, and changed pinned defaults.
- Documentation examples and installation checks pass.
- Implementation is merged to `master`, and its worktree is removed.

Completion evidence will be recorded in `docs/VALIDATION.md` and the tracking issue.
