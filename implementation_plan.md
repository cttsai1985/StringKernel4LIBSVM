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
