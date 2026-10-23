# Repository instructions

- Read `docs/DESIGN.md`, `docs/COMPATIBILITY.md`, and `docs/DEVELOPMENT.md` before
  implementation. Check existing branches and worktrees for duplicate work.
- Use one worktree per open PR. Reuse an appropriate worktree before creating
  another. Keep scratch output under ignored build directories or `/tmp`.
- Use TDD for behavior changes. Observe an independently specified regression
  test fail, implement the fix, then run it successfully before refactoring.
- Preserve the CLI and mathematical contracts in `docs/COMPATIBILITY.md` and
  `docs/KERNELS.md`. Add pin tests for production defaults and dependency changes.
- Use values, standard containers, and explicit configuration. Library code
  reads no clock. Pass any future timestamps as arguments.
- Keep vendored LIBSVM files immutable against their upstream pin. Update source,
  provenance, independent digest tests, and notices together when upgrading.
- Keep test dependencies behind `BUILD_TESTING`. Production builds must work
  without Python or network access.
- Keep `svmsk_analyser.cpp` unchanged and excluded from supported targets unless
  a separate design restores its missing dependencies and public API integration.
- Run the applicable commands in `docs/DEVELOPMENT.md`. Verify real matrices and
  predictions before claiming success. Label measurements `MEASURED` and estimates
  `ESTIMATED`. Do not infer numerical correctness from compilation alone.
- Finish authorized work through issue, PR, current-commit CI, and merge to
  `master`. Close the issue with evidence and remove the worktree after merge.
  Never use `git reset --hard` or `git stash` in the operator's working tree.
- Retain upstream notices. Do not assign the historical project a new license.
