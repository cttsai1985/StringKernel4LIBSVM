# Contribute a change

Read [the design](docs/DESIGN.md) and [compatibility contract](docs/COMPATIBILITY.md)
before changing behavior. Use [development instructions](docs/DEVELOPMENT.md) for
the build, test, formatting, sanitizer, and coverage commands.

## Fix behavior with a regression test

1. Open or update an issue with the input, expected result, and acceptance checks.
2. Check existing branches and worktrees for the same work. Use one worktree per PR.
3. Add a test that exposes the defect. Run it and record its expected failure.
4. Make the smallest change that passes the test.
5. Refactor with the tests passing. Check adjacent code for the same defect.
6. Run the full development suite and the checks affected by the change.
7. Update the documentation and link the evidence in the PR.

Numerical tests need an independent expected result. Use hand calculations or a
separate feature-map construction rather than another copy of the production
loop. Include cases that reject wrong degree weights, position handling, or
normalization. A successful process exit alone does not prove a correct matrix or
prediction. Integration tests must inspect the actual output.

Pin production defaults and dependency versions in tests. Label reported
measurements `MEASURED` and estimates `ESTIMATED`. Give the command, environment,
input, and comparison for a performance claim. Correctness fixtures do not
establish biological prediction quality.

## Keep dependencies reproducible

Keep `third_party/libsvm` identical to its recorded upstream commit. To update
LIBSVM, review the upstream change, replace the complete pinned source set, and
update `UPSTREAM.json` and the independent expectations in
`tests/test_dependency_pins.py`. Retain the upstream copyright notice.

Catch2 is a test-only dependency. Update its version, archive digest, independent
pin test, and documentation together. Preserve production configuration with
`BUILD_TESTING=OFF` so it needs neither Python nor network access.

## Finish the change

The PR must pass CI on its current commit before merge to `master`. Close the
tracking issue with evidence after merge, then remove the implementation
worktree. Keep a regression gate that fails if the fixed behavior returns.
Do not rewrite an operator's working tree with `git reset --hard` or `git stash`.

The owner has not selected a project license. Do not add or imply a new license
as part of an unrelated contribution.
