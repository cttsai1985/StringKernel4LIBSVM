# LIBSVM bootstrap cleanup plan

Status: implementation and local/fresh-archive verification complete; PR #4 targets
`dev` for independent root review and exact-head hosted CI. Master promotion and
issue #3 closeout remain with the root. Preserve ordered commits. Previous modernization evidence remains in `docs/VALIDATION.md` and Git.

## Goal and risk

Remove the committed upstream LIBSVM copy and fetch the same unmodified official
3.37 commit into an ignored directory through a small standard-library script.
Risk: moderate, because dependency availability and build setup change; solver,
CLI and mathematical output stay compatible; sample data and historical analyzer
move to clear directories without byte changes.

## Scope (exact blast radius)

Modify: `implementation_plan.md`, `AGENTS.md`, `.gitignore`, `.gitattributes`,
`CMakeLists.txt`, `tests/CMakeLists.txt`, `tests/test_dependency_pins.py`,
`.github/workflows/ci.yml`, `README.md`, `CONTRIBUTING.md`,
`THIRD_PARTY_NOTICES.md`, `docs/DESIGN.md`, `docs/DEVELOPMENT.md`,
`docs/VALIDATION.md`, `docs/COMPATIBILITY.md`, `sample.bat`,
`tests/test_libsvm_integration.py`, `tests/test_options.cpp`,
`StringKernelConverter.cpp`.

Add: `cmake/libsvm.json`, `cmake/Libsvm.cmake`, `scripts/fetch_libsvm.py`,
`tests/test_fetch_libsvm.py`.

Move unchanged data: `data.txt` -> `examples/regression/train_sequences.tsv`,
`test_data.txt` -> `examples/regression/test_sequences.tsv`.
Move historical bytes unchanged: `svmsk_analyser.cpp` ->
`legacy/paired_sequence_analyzer.cpp`.
Rename internal parser: `GetLongOpt.h` -> `CliOptions.h`,
`GetLongOpt.cpp` -> `CliOptions.cpp`; update active includes and coverage paths.

Remove tracked copies: `third_party/libsvm/svm.cpp`, `third_party/libsvm/svm.h`,
`third_party/libsvm/svm-train.c`, `third_party/libsvm/svm-predict.c`,
`third_party/libsvm/README`, `third_party/libsvm/COPYRIGHT`,
`third_party/libsvm/UPSTREAM.json`.

## Changes and architecture

1. Keep a small tracked manifest at `cmake/libsvm.json` with exact commit
   `6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb`, official URLs, and existing SHA-256
   digests. Retain independent expected pins in tests.
2. First write failing behavioral tests for a standard-library bootstrap:
   successful verified fetch, offline reuse, incomplete download recovery,
   corruption rejection, bad download hash, and network failure. Mock only the
   external HTTP boundary; exercise real bytes, digests, filesystem, and CLI.
3. Implement `python3 scripts/fetch_libsvm.py` with optional `--source-dir`.
   Download missing files into temporary staging, verify all required bytes
   before publishing files, reuse a fully verified cache without network access,
   and report actionable errors. Keep fetched bytes ignored at
   `third_party/libsvm`; do not create a private fork.
4. Add a pure CMake preflight which validates required files and their SHA-256
   hashes against the manifest and reports the bootstrap command when needed.
   Configuration/build after bootstrap must need neither Python nor network.
5. Bootstrap before configure in CI. Document the setup command first on Unix
   and Windows, explain offline reuse, and preserve installation of the fetched
   COPYRIGHT beside the existing notice. Update contributor rules accordingly.

6. Relocate examples and historical source to the scoped paths; pin their original
   byte digests in a regression test. Update sample workflow, active parser names,
   docs, formatting selection, and integration paths without altering numerics.
7. Verify a fresh Git archive: no tracked LIBSVM bytes, actionable failure before
   bootstrap, and successful pinned fetch/build/test afterward.

## Risks and alternatives

Initial setup needs Python 3.9+ and official-source network access. Verified local
sources can be carried to an offline machine. Existing corrupt files fail rather
than being silently trusted. A failed fetch must leave no new usable partial
cohort. No LIBSVM fork is needed because the solver has no local patches.

## Verification and completion

- Record RED then GREEN fetcher tests and independent dependency pins.
- Verify a fresh production configure fails with an actionable bootstrap error.
- Bootstrap official bytes, verify the complete development suite, real solver
  matrices/predictions, and installed artifacts.
- Verify production with testing/Python discovery/network fetching disabled.
- Run sanitizers, formatting, and project coverage (98% line floor).
- Obtain current-head Linux/macOS/Windows/sanitizer/coverage/format hosted CI and
  independent root review; open a PR to `dev` linking the issue without closing
  it until promotion to `master`. Preserve ordered commits without squashing.
- After authorized master promotion, close the issue with evidence and remove
  the worktree. The implementation agent must stop short of merging for review.


---

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

Spectrum uses the frequency-count feature map of contiguous, position-independent k-mers; its raw kernel is the dot product of those counts. SSK uses all ordered length-k subsequences with gap weighting. The implementation may use the cosine-equivalent gap-only feature scaling `lambda^(span-k)` to avoid a common underflow factor; documentation must state that scaling and compare normalized values to the paper's `lambda^span` definition. For either new fixed-length kernel, `k > sequence length` gives an empty feature vector and follows the existing zero-self normalization convention. WDSK still uses all available substring orders up to the sequence length, and DSK does not use `k` in its score. Precompute/check all sequence self-values before writing; reject non-finite results without partial matrix output. Never silently turn a nonempty SSK feature vector into a zero diagonal.

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

---

# Kernel separation analysis

Status: locally implemented and verified; awaiting current-head hosted review and ordered promotion.

## Goal

Add a reproducible, approachable analysis showing how positional, weighted-degree, spectrum, and subsequence kernels separate three small labeled sequence datasets with distinct alignment, movement, and gap patterns.

## Exact blast radius

- Modify `implementation_plan.md` — this plan and measured evidence.
- Add `examples/analysis/aligned-motif.tsv` — tiny equal-length aligned motif set.
- Add `examples/analysis/moved-motif.tsv` — same motif at different positions.
- Add `examples/analysis/gapped-motif.tsv` — ordered motif with insertions/gaps.
- Add `scripts/analyze_separation.py` — stdlib-only reproducible runner; invoke the built converter and actual LIBSVM tools, independently calculate class-pair similarity summaries, and emit a self-contained SVG/Markdown report to a caller-selected output directory.
- Add `tests/test_analysis.py` — check generated report content, deterministic summaries, output safety, and actual solver predictions using supplied executable paths.
- Modify `tests/CMakeLists.txt` — register analysis behavior/integration test with built executable paths.
- Add `docs/ANALYSIS.md` — explain dataset construction, similarity heatmaps, results, limitations, and exact rerun command.
- Add `docs/analysis/similarities.svg` — published deterministic figure; the integration test pins the generated SVG bytes to this artifact.
- Modify `README.md` — link the analysis beside kernel-selection guidance.

## Design and sequence

1. Define the three miniature datasets with transparent labels and held-out examples; use no external data, package, network, statistical claim, or newly inferred biological alphabet.
2. Write the test first. It invokes the runner against supplied converter/trainer/predictor paths and asserts expected labels, independently recomputed pairwise similarities, deterministic report shape, SVG validity, and no writes outside the selected output directory. Record the initial RED result.
3. Implement a standard-library analysis runner. Converter output remains the source of actual solver Gram matrices; a separate tiny feature-map oracle calculates per-class within/between similarity summaries. Train/predict through the real LIBSVM binaries on designated supported-kernel/dataset combinations. Fail with useful diagnostics on process or artifact errors.
4. Generate a readable SVG heatmap and Markdown summary using explicit kernel and parameter settings. Do not add plots, ML, or other dependencies; commit only the SVG intentionally published as documentation, not build outputs.
5. Document what the examples demonstrate and what they cannot establish. Distinguish observations on these constructed records from accuracy, generalization, or comparative speed claims.
6. Run focused RED/GREEN analysis checks, then the repository's development/sanitizer/coverage/format checks and exact-head hosted CI. Obtain independent review; merge feature PR to `dev`, promote by a separate PR to `master`, close the issue with evidence, and remove both worktrees.

## Risks and trade-offs

Synthetic records make the kernel differences legible but are not representative datasets. Results must be framed as geometry demonstrations and correctness examples. Keeping all dataset tokens tiny bounds the explicit subsequence oracle to test-only inputs. SVG generation avoids a plotting dependency but needs a test that validates key data labels and heatmap dimensions. The runner writes only beneath an explicitly selected output directory and does not overwrite source data.

## Acceptance and proof

- Three datasets each expose the intended alignment/movement/gap distinction.
- Actual converter, LIBSVM train, and predict run successfully; predictions match the held-out labels selected by the fixture design.
- Similarity numbers match an independently specified feature-map calculation; within-class/between-class summaries are labeled as this fixture's observations.
- The SVG and Markdown are regenerated deterministically and tested for content.
- No runtime dependencies, kernel API changes, production behavior changes, public data provenance assumptions, accuracy claims, or speed claims are added.
- All current-head hosted gates pass; promotion is merged to master and issue closed with artifact-backed evidence.

RED evidence: after adding the end-to-end report test and CTest registration, direct execution failed because `scripts/analyze_separation.py` did not yet exist (expected missing implementation). No production code was added before this check.


## Local completion evidence

- RED: the first analysis integration test failed because `scripts/analyze_separation.py` did not yet exist.
- MEASURED: development, ASan/UBSan, and coverage configurations each passed all 50 CTest checks after implementation.
- MEASURED: project-only line coverage remained 391/398 (98.2%); the existing 98% floor passed.
- MEASURED: all three actual converter matrices matched the separate feature-map oracle; real LIBSVM training and prediction recovered `1, -1` on each held-out pair.
- MEASURED fixture means (within-class, between-class): aligned weighted degree 0.700/0.425; moved spectrum 0.500/0.000; gapped subsequence 0.013/0.000. These describe only the constructed four-record training sets.
- MEASURED: regenerated SVG SHA-256 equals the published `docs/analysis/similarities.svg` digest `e22ff26b944e4fbcdab72dcac911406f3194678e567dbc6b4004fcd6f5bb90c1`.
- Python syntax compilation and `git diff --check` passed. Hosted CI and independent review remain required before promotion.
