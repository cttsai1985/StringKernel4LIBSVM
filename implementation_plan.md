# Spectrum and subsequence kernels

Status: approved scope; implementation is proceeding in ordered test-first commits.

## Goal

Add two common, normalized string kernels to the supported C++17 library and CLI: a fixed-length spectrum kernel and a fixed-length decay-weighted subsequence kernel (SSK). Preserve existing positional and weighted-degree contracts and defaults. Improve the README with a small runnable LIBSVM classification example and an approachable kernel-selection guide.

## Exact blast radius (17 paths)

1. `implementation_plan.md` — this scoped plan and commit/verification record.
2. `docs/DESIGN.md` — update API and CLI design contracts.
3. `StringKernels.h` — explicit enum IDs, decay configuration, public helper declarations.
4. `StringKernels.cpp` — raw/self/normalized spectrum and SSK calculations.
5. `CliOptions.cpp` — named options and parameter validation, preserving `-t 0/1`.
6. `LibsvmWriter.cpp` — generic self-kernel precomputation and new dispatch.
7. `tests/test_kernels.cpp` — independent mathematical oracles and edge cases.
8. `tests/test_options.cpp` — legacy pins and new CLI behavior.
9. `tests/test_writer.cpp` — exact matrices and pre-output failures.
10. `tests/test_cli.py` — process-level CLI coverage.
11. `tests/test_libsvm_integration.py` — both kernels through real solver training/prediction.
12. `examples/string-kernels/spectrum-train.tsv` — tiny labeled fixture.
13. `examples/string-kernels/spectrum-test.tsv` — held-out fixture.
14. `docs/KERNELS.md` — formulas, use cases, and analytical costs.
15. `docs/COMPATIBILITY.md` — stable old values and new options.
16. `README.md` — newcomer guide and runnable examples.
17. `docs/VALIDATION.md` — RED/GREEN and measured completion evidence.

No build workflow, CMake target, dependency, license, or legacy analyzer change is expected. If existing targets fail to discover tests or fixtures, stop and get the blast radius revised before editing more paths.

## Stable contracts and mathematics

Keep `distance=0`, `weighted_degree=1`, default weighted degree, degree 5, existing byte/case semantics, unequal-length behavior, and LIBSVM output stable. Add `spectrum=2` and `subsequence=3`; support readable `--kernel positional|weighted-degree|spectrum|subsequence`, retaining existing numeric flags. `-k` is positive and means max weighted-degree order, exact spectrum k-mer length, or exact SSK subsequence length. `--decay` is valid only for SSK and must be finite in `(0,1]`, default `0.5`.

Spectrum uses the frequency-count feature map of contiguous, position-independent k-mers; its raw kernel is the dot product of those counts. SSK uses all ordered length-k subsequences with gap weighting. The implementation may use the cosine-equivalent gap-only feature scaling `lambda^(span-k)` to avoid a common underflow factor; documentation must state that scaling and compare normalized values to the paper's `lambda^span` definition. `k > sequence length` gives an empty feature vector and follows the existing zero-self normalization convention. Precompute/check all sequence self-values before writing; reject non-finite results without partial matrix output. Never silently turn a nonempty SSK feature vector into a zero diagonal.

## Ordered commits and verification

1. Commit this plan and open the issue-linked PR to `dev` before code.
2. Add independently specified behavior tests, including pins for old IDs/defaults and short-string spectrum/SSK feature-enumeration oracles. Run them and record RED.
3. Implement spectrum; run focused and full kernel checks.
4. Implement SSK with explicit numeric/error behavior; rerun independent oracles and boundary cases (repetition, gaps, tiny positive decay, decay 1, and `k > length`).
5. Add CLI/writer support and exact normalized matrix tests, including no output before invalid/non-finite checks.
6. Exercise both kernels with real LIBSVM using the documented fixture and independent matrix expectations; expected labels are correctness checks, not quality claims.
7. Update design, compatibility, kernel and README docs after behaviors are verified. Record only measurements actually collected.
8. Run the existing `docs/DEVELOPMENT.md` suite: build/CTest, CLI and solver integration, format, sanitizers, dependency pin checks, install workflow, coverage floor, and Linux/macOS/Windows hosted CI. Independent review must pass on the exact PR head before master promotion.

## Examples and limits

The fixture uses `ACG`/`ACT` as positive and `GGA`/`GGT` as negative training examples, with `ACC` and `GGC` held out. Spectrum `k=2` recognizes shared `AC` and `GG` features. For SSK, `ACG` and `AG` share the ordered gapped subsequence `AG`; contiguous spectrum `k=2` does not. Run both through the actual converter and LIBSVM train/predict harness. State that this is a correctness fixture, not a generalization benchmark.

Explain that DSK compares matching positions, WDSK compares same-position substrings across orders, spectrum compares exact contiguous k-mers regardless of position, and SSK compares ordered gapped patterns while discounting gaps. Inputs remain whitespace-delimited byte tokens, not FASTA; no alphabet validation, alignment, tuning, or cross-validation is implied.

Derived complexity: with mean token length `L`, SSK pair evaluation is `O(k L^2)` and a dense `N x N` matrix is `O(N^2 k L^2)` before solver cost; mark this analytical bound, not a benchmark. Spectrum is expected `O(k(Lx+Ly))` per pair with count maps; label this as an algorithmic/implementation estimate, not measured timing. Do not make speed or accuracy claims without a benchmark that can lose.
