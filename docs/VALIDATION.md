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
