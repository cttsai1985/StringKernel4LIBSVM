# Kernel and matrix reference

## Weighted-degree string kernel

For requested degree `d > 0`, the weight of substring order `j` is

```text
beta(j) = 2 * (d - j + 1) / (d * (d + 1)),  1 <= j <= d
```

Let `L` be the shorter sequence length. The raw kernel is

```text
raw(a,b) = sum over j = 1..min(d,L) of
           beta(j) * sum over i = 0..L-j of [a[i:i+j] == b[i:i+j]]
```

The indicator contributes one for an exact match at the same starting position.
The normalized result is

```text
K(a,b) = raw(a,b) / sqrt(raw(a,a) * raw(b,b))
```

The implementation computes each self-kernel from the sequence's entire length.
It returns zero when either self-kernel is zero. Dataset input rejects empty
sequences, but the lower-level weighted-degree functions accept empty views.
The normalization helper rejects negative or nonfinite inputs.

The requested degree controls the weights even when a sequence is shorter than
`d`. Increasing `d` does not add nonexistent substrings or allocate storage
proportional to `d`. The default degree is `5`.

At degree `2`, `beta(1) = 2/3` and `beta(2) = 1/3`. Therefore
`raw(AC,AC) = 5/3`, `raw(AC,AD) = 2/3`, and `K(AC,AD) = 2/5`.
These fractions are mathematical specifications, not benchmark measurements.

### Position and unequal lengths

This kernel compares positions. A matching substring at another position does
not contribute to that position's score. No alignment, case conversion, alphabet
validation, or Unicode normalization occurs.

For unequal lengths, this project preserves its historical extension: compare
substrings within the shared positions and normalize with both complete
self-kernels. Sequence suffixes outside the shared span contribute only to their
sequence's self-kernel. This extension is part of the project's compatibility
contract.

The aligned weighted-degree kernel is described in section 4.2.3 of
[Sonnenburg et al., Large Scale Multiple Kernel Learning](https://www.mlsec.org/docs/2007-lskm.pdf).

## Spectrum kernel, converter type `2`

For a positive `k`, each coordinate is a token of `k` consecutive bytes and its
feature value is the number of times that token occurs:

```text
phi_u(x) = count of contiguous k-byte substrings equal to u
raw(x,y) = sum over u of phi_u(x) * phi_u(y)
K(x,y) = raw(x,y) / sqrt(raw(x,x) * raw(y,y))
```

This is a feature-map dot product, so it is positive semidefinite. Occurrences
are counted with multiplicity; their positions do not matter. For example, `ACG`
and `GAC` share the `AC` 2-mer even though it moved, while `ACG` and `AG` have
no contiguous 2-mer in common. If `k` exceeds a sequence length its feature map
is empty, and the normalized result is zero by the zero-self convention.

The `--kernel spectrum -k K` selector uses a positive exact k-mer size. The
numeric alias `-t 2` is stable. An independent implementation of the spectrum
kernel is available in the [KeBABS sequence-kernel package](https://www.bioconductor.org/packages/release/bioc/html/kebabs.html).
This project implements the published count-map definition directly and does
not depend on or copy that package.

## String subsequence kernel, converter type `3`

For each byte string `u` of length `k`, the paper's feature value sums over all
strictly increasing embeddings `I` of `u` into `x`, weighting each embedding by
`lambda^(span(I))`, where `span(I) = last(I) - first(I) + 1`. The SSK is the dot
product of those feature vectors, then cosine normalized. Its feature map makes
it positive semidefinite. A smaller `lambda` reduces the contribution of
embeddings with gaps. Occurrences with the same symbols add together, including
repeated symbols.

The implementation uses the cosine-equivalent gap-only feature weight
`lambda^(span(I)-k)`. This divides every feature coordinate from the paper's
definition by the same positive `lambda^k`, so the normalized kernel is
unchanged while avoiding a common `lambda^(2k)` raw-kernel factor. It does not
change which embeddings contribute. `--kernel subsequence -k K --decay LAMBDA`
selects this kernel; `-t 3` is the stable numeric alias. `LAMBDA` must be finite
and in `(0, 1]`; the default is `0.5`.

For `ACG` and `AG` with `k=2` and `lambda=0.5`, the only shared feature is `AG`.
Its gap-only feature value is `0.5` in `ACG` and `1` in `AG`. The `ACG` self
value is `1^2 + 0.5^2 + 1^2 = 2.25`, so the normalized kernel is `1/3`. For
spectrum `k=2`, `AG` has an empty feature map and the result is zero.

The dynamic program adds one matched position pair at each order and discounts
the gaps between consecutive matches. With sequence lengths `n` and `m`, it uses
`O(k*n*m)` operations and `O(n*m)` working memory for a pair. A dense `N` by `M`
precomputed matrix therefore has analytical SSK work `O(N*M*k*n*m)`, plus
per-sequence self-kernel precomputation. For equal length `L`, this is
`O(N*M*k*L^2)`. These are derived operation bounds, not timing measurements.
SSK can be substantially more costly than exact spectrum counts; use it when
gapped ordered matches fit the problem and sequence lengths are manageable.

The writer computes all sequence self-kernels before writing a matrix row. If
the SSK dynamic program or its final value exceeds the supported numeric range,
conversion fails instead of emitting a non-finite value or silently replacing a
nonempty diagonal with zero. `k > length` remains a valid empty-feature case.

For both new kernels, an empty feature map returns normalized zero. Identical
sequences with a nonempty feature map normalize to one. Inputs are the existing
case-sensitive byte tokens: no DNA/protein alphabet, FASTA parser, or alignment
is inferred.

## Choosing a kernel

| Kernel | What matches | A useful fit |
| --- | --- | --- |
| DSK (`0`, `positional`) | Same byte at the same position | Already aligned, equal-length strings where position matters. |
| WDSK (`1`, `weighted-degree`) | Same-position substrings across orders up to `k` | Aligned strings where short and longer local matches matter. |
| Spectrum (`2`) | Exact contiguous k-mers anywhere | Position can shift, but each k-mer must match exactly. |
| SSK (`3`) | Ordered k-byte subsequences, including gaps | Small local insertions/deletions should preserve partial similarity; account for its higher cost. |

For a small readable example and matching train/predict commands, see the
[README](../README.md#choose-a-kernel-and-run-a-small-classification-example).

Primary definitions: [Leslie, Eskin & Noble, The spectrum kernel (PSB 2002)](https://cseweb.ucsd.edu/~eeskin/papers/spectrum-psb02.pdf) and
[Lodhi et al., Text classification using string kernels (JMLR 2002)](https://www.jmlr.org/papers/volume2/lodhi02a/lodhi02a.pdf).
For additional reference implementations, see KeBABS's [spectrum kernel package](https://www.bioconductor.org/packages/release/bioc/html/kebabs.html)
and Shogun's [subsequence string-kernel source](https://github.com/shogun-toolbox/shogun/blob/8f01b2b9e4de46a38bf70cdb603db75ebfd4b58b/src/shogun/kernel/string/SubsequenceStringKernel.cpp).
These links are references only; this project contains its own feature-map
implementation and adds no runtime dependency on either package.

## Positional similarity, converter type `0`

The legacy name `DSK` denotes the fraction of matching positions:

```text
DSK(a,b) = number of positions i with a[i] == b[i] / sequence length
```

All reference and target sequences must have the same nonzero length. The writer
validates that condition before emitting a row. The kernel is positional
similarity and does not compute edit distance. Its result does not depend on
degree, but configuration still requires a positive degree.

## Precomputed matrix format

Each output row has the form

```text
label 0:row_serial 1:K(target,reference_1) ... N:K(target,reference_N)
```

`N` is the number of reference records. All columns appear, including explicit
zeros. Column order is the original reference order. Serial numbers start at
one in each output file. LIBSVM uses training serials to identify reference rows
and ignores prediction serials.

Training conversion uses the training dataset as both reference and target,
producing a square matrix. Prediction conversion uses the original training
dataset as reference and the test dataset as target. Its row count equals the
number of test records, while its column count remains `N`.

Both conversions must use the same kernel configuration and reference order.
Reordering the test targets changes row order only. Reordering the reference
changes the meaning of every column and requires a newly trained model.

Values use the classic locale and round-trip double precision. Labels remain
finite doubles, including regression labels. Validation errors occur before
matrix output. An I/O error during output can leave a partial matrix and returns
a failing CLI status.

The converter already normalizes its kernels. Applying `svm-scale` to these
matrices as ordinary features changes the kernel and can corrupt the serial
column. LIBSVM's
[precomputed-kernel specification](https://github.com/cjlin1/libsvm/blob/v337/README#precomputed-kernels)
describes the solver's input contract.
