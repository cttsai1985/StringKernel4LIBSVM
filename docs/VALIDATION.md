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
