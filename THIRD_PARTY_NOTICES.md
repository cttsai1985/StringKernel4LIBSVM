# Third-party notices

## LIBSVM

The build includes unmodified source from
[LIBSVM 3.37, tag `v337`](https://github.com/cjlin1/libsvm/tree/v337), commit
`6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb`.

Copyright (c) 2000-2023 Chih-Chung Chang and Chih-Jen Lin. All rights reserved.

The full redistribution conditions and disclaimer are in
[the pinned upstream COPYRIGHT](https://github.com/cjlin1/libsvm/blob/6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb/COPYRIGHT).
The dependency bootstrap fetches that complete file into ignored
`third_party/libsvm/COPYRIGHT`; installation copies it beside this notice.
The exact files, source URLs, and SHA-256 digests are recorded in
[`cmake/libsvm.json`](cmake/libsvm.json).
In an installation, the complete `COPYRIGHT` file is beside this notice.

## Catch2

Tests use [Catch2 v3.16.1](https://github.com/catchorg/Catch2/tree/v3.16.1), licensed
under [BSL-1.0](https://github.com/catchorg/Catch2/blob/v3.16.1/LICENSE.txt).
Catch2 is downloaded only when tests are enabled and no local source override is
provided. It is not linked into the converter or installed project library.

## Project license

The historical StringKernel4LIBSVM repository does not declare a project license.
This modernization retains third-party notices and does not assign a new license
to project-owned code.
