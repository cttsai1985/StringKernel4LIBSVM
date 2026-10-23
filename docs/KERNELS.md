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
