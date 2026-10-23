# Third-party notices

## LIBSVM

The build includes unmodified source from
[LIBSVM 3.37, tag `v337`](https://github.com/cjlin1/libsvm/tree/v337), commit
`6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb`.

Copyright (c) 2000-2023 Chih-Chung Chang and Chih-Jen Lin. All rights reserved.

The full redistribution conditions and disclaimer are in
[`third_party/libsvm/COPYRIGHT`](third_party/libsvm/COPYRIGHT). This file remains
part of the source distribution and is installed with this notice. The vendored
files and their source URLs and SHA-256 digests are recorded in
[`third_party/libsvm/UPSTREAM.json`](third_party/libsvm/UPSTREAM.json).

## Catch2

Tests use [Catch2 v3.16.1](https://github.com/catchorg/Catch2/tree/v3.16.1), licensed
under [BSL-1.0](https://github.com/catchorg/Catch2/blob/v3.16.1/LICENSE.txt).
Catch2 is downloaded only when tests are enabled and no local source override is
provided. It is not linked into the converter or installed project library.

## Project license

The historical StringKernel4LIBSVM repository does not declare a project license.
This modernization retains third-party notices and does not assign a new license
to project-owned code.
