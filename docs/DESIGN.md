# Design: a tested sequence-kernel converter with official LIBSVM

## Supported workflow and boundary

The supported product accepts whitespace-delimited sequence/label records,
produces a precomputed matrix, and ships official LIBSVM training and prediction
programs. A project-owned C++17 library handles sequence mathematics and data
formatting. The solver remains unmodified. This boundary preserves the original
project purpose and avoids a private solver fork.

`StringKernelConverter -train train.txt [-test target.txt] [-k degree] [-t type]`
continues to work. Default type is weighted degree (`1`), default degree is `5`;
type `0` is equal-length positional similarity. Add `--help` and `--version`.
Missing/unknown/duplicate options and incomplete or out-of-range integer tokens
are errors. Negative degree values are parsed as values and rejected, not
mistaken for option names. Help/version are standalone commands.

`CMakeLists.txt` owns project version `1.0.0`; the executable reports that version
and the bundled LIBSVM version, and tests pin both. This is the modernization's
initial explicit project version, not a claim of historical release continuity.

## Data model and API

All project-owned declarations live in `namespace stringkernel`. Public headers
are self-contained, use scoped names, and contain no global `using namespace`.

```cpp
enum class KernelType { distance = 0, weighted_degree = 1 };
struct KernelConfig {
    KernelType type = KernelType::weighted_degree;
    int degree = 5;
    void validate() const;
};
double degree_weight(int degree, std::size_t order);
double weighted_degree(std::string_view a, std::string_view b, int degree);
double weighted_degree_self(std::string_view sequence, int degree);
double normalize_kernel(double value, double self_a, double self_b);
double kernel_value(std::string_view a, std::string_view b,
                    const KernelConfig& config);

struct Sample { std::string sequence; double label; };
using Dataset = std::vector<Sample>;
Dataset read_dataset(std::istream& input, std::string_view source = "<stream>");
Dataset load_dataset(const std::filesystem::path& path);
void write_precomputed(std::ostream& output, const Dataset& reference,
                       const Dataset& target, const KernelConfig& config = {});

struct Options {
    std::filesystem::path train;
    std::optional<std::filesystem::path> test;
    KernelConfig kernel;
    bool help = false;
    bool version = false;
};
Options parse_options(const std::vector<std::string_view>& arguments);
std::string_view usage_text();
```

The argument vector excludes the program name. The library never reads the clock
or changes process locale. Ownership is entirely in values/standard containers;
there are no writable global datasets or normalization caches. Kernel parameters
are explicit at every calculation. CLI parsing is internal to the executable;
the installed library headers are kernels, datasets, and the writer.

## Input and failure contracts

Each nonempty record contains exactly a sequence token and a finite double label.
Both spaces and tabs delimit fields; whitespace-only records are ignored. Support
LF, CRLF, and a final record without newline. Parse numbers in the classic locale,
consume the entire numeric token, and reject trailing fields. Sequences are
case-sensitive nonwhitespace byte strings, without a newly imposed amino-acid
alphabet, implicit uppercase conversion, or automatic alignment.

Use `getline`, with no fixed record-size buffer. Reject empty datasets. A malformed
nonempty record raises a diagnostic identifying source and one-based line number.
File-open/read failures include their path. Validate datasets and configuration
before writing the first row; direct library callers receive the same validation
as file callers, including finite labels and nonempty sequence values.

Library contract violations use `std::invalid_argument`; malformed/unreadable
input and output failures use `std::runtime_error`. Main catches standard
exceptions, writes the diagnostic to stderr, and returns failure. Usage/version
and successful conversion return success. Conversion stdout is matrix data only.
Flush/check output at completion so redirected write failures cannot report
success. A late I/O failure can leave partial output and must return failure.

## Mathematics and compatibility

For requested degree `d`, the order weight is

`beta(j) = 2 * (d - j + 1) / (d * (d + 1)), 1 <= j <= d`.

The raw weighted-degree kernel sums `beta(j)` for every same-position exact
substring match of order `j`. The normalized kernel is

`K(a,b) = raw(a,b) / sqrt(raw(a,a) * raw(b,b))`.

Zero self-kernels yield zero rather than NaN. Reject nonfinite or negative inputs
to the normalization helper. Avoid multiplying self-kernels before taking square
roots, which can overflow even when the normalized result is finite.

As an analytical fixture, at degree `2`, `raw(AC,AC) = 5/3`,
`raw(AC,AD) = 2/3`, and the normalized `K(AC,AD) = 2/5`.
These are exact mathematical specifications, not measured benchmark results.

Preserve the legacy unequal-length WDSK extension: compare only shared positions,
but use each entire sequence's self-kernel. Preserve the requested degree in
weights when degree exceeds input length. Accept positive degrees representable
as `int`, using floating-point denominators to avoid signed integer overflow;
never allocate a vector proportional to an arbitrary requested degree.

A bounded matching-prefix loop avoids substring allocation: for each starting
position, extend matching prefixes only through `min(d, remaining positions)` and
add their degree weights. Its worst-case operation bound is analytically
`O(L * min(d,L))` per pair, where `L` is the shorter length. This is an algorithmic
bound, not a measured speedup. Self-kernels can sum `(length-j+1)*beta(j)` directly.

DSK is positional similarity, not edit distance. For equal nonempty lengths,
return matching positions divided by sequence length. Reject unequal lengths
across all reference and target records before matrix output; this prevents the
legacy reference-only normalization from producing asymmetric training matrices.
The dormant RBF-over-similarity code is not part of the supported API.

Use round-trip double precision and classic-locale decimal points. This preserves
mathematical values, not the old byte format with six decimal places. Integer
class labels and finite regression labels are both retained without truncation.

## Matrix generation and resource behavior

Reference dataset order defines columns permanently. Train/train output is square;
target/train output has a row per target and a column per reference. Every row is
`label 0:row_serial 1:K(target,ref1) ... N:K(target,refN)`, including explicit zeros.
Serials are one-based; prediction row serials may exceed reference size because
LIBSVM ignores prediction serials. Train serials are in the reference range.

Validate all records first. For WDSK, precompute self values once per dataset per
call; when reference and target are the same object, reuse those values. Stream
rows into the supplied output stream instead of storing the matrix. Caller-owned
input storage is linear in sequence data; auxiliary self storage is linear in
record count, independent of the number of matrix cells. Fail before conversion
if reference count exceeds the integer-index limit supported by LIBSVM.

## Build, dependencies, and installation

Require CMake `3.25` or newer. Presets use a platform-native generator and explicit
build/test configuration, supporting both single- and multiple-configuration
generators. Keep project warning settings and sanitizer/coverage instrumentation
away from upstream-source warning policy. Build upstream CLI `.c` files as C++,
matching the official Makefile. Offer sanitizer and coverage configurations only
on supported compiler/platform combinations and fail clearly on invalid requests.

Fetch LIBSVM tag `v337`, exact commit
`6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb`, with upstream sources, README, and
COPYRIGHT unchanged. `cmake/libsvm.json` records the version, commit, official
source URLs, and SHA-256 for every required file. The standard-library
`scripts/fetch_libsvm.py` bootstraps ignored `third_party/libsvm`, verifies cached
bytes before offline reuse, and stages missing downloads until the entire cohort
is verified. Modified cached bytes and bad downloads fail with actionable errors.
No downloaded upstream source is committed and no solver fork is maintained.

CMake independently hashes every required cached file on configuration. Missing
or modified sources fail with the bootstrap command; CMake never runs Python or
uses the network to prepare LIBSVM. Tests independently pin the expected commit
and hashes, so editing source and manifest together is detected. Python 3.9+ is
required for bootstrap and tests, while production builds after bootstrap are
independent of Python and network access.

Catch2 `v3.16.1` is test-only and fetched by archive digest. An offline source
override is documented. Python standard-library test drivers, Python discovery,
and Catch2 setup are all gated behind `BUILD_TESTING`. After bootstrap, production
configuration with `-DBUILD_TESTING=OFF` needs neither Python nor network access.

Install the converter, official trainer/predictor, project library/public headers,
and upstream notices. An integration test installs into a temporary prefix and
runs that installation. Keep project license choice explicit; this work does not
silently assign a license to the historical project.

## Tests and delivery gates

The owner requested pstack and PUA as workflow guidance. Apply pstack's
[Test Behavior, Not Implementation](https://github.com/cursor/plugins/blob/main/pstack/skills/principle-test-behavior-not-implementation/SKILL.md)
by testing independent mathematical expectations and the public CLI.
Apply [Prove It Works](https://github.com/cursor/plugins/blob/main/pstack/skills/principle-prove-it-works/SKILL.md)
by running real solver and installed-artifact workflows, and its
[shipping verification](https://github.com/cursor/plugins/blob/main/pstack/skills/poteto-mode/playbooks/shipping.md)
by obtaining an independent verdict on the current patch before merging. Use
available Codex collaborators and CLI tools for these checks; no Cursor cloud
agent or unavailable harness-specific tool is assumed. PUA supplies persistence,
root-cause investigation, alternate approaches after repeated failures, and
evidence-backed issue closeout. Structural refactoring and intentional behavior
corrections remain distinguished in the compatibility notes and test evidence.

Numerical tests build a separate feature map keyed by `(position, order, token)`
and compare weighted dot products to production results. Include hand calculations,
bounded exhaustive small alphabets, positional shifts, repeated symbols, degree
boundaries, unequal lengths, symmetry, diagonal, and zero/invalid normalization.
This oracle must be able to reject a wrong weighting or normalization formula.

Writer tests assert exact dimensions/column order/zeros, repeated calls with changed
configuration, dataset permutations, single-target equivalence, input rejection
before output, classic-locale output, precision, and failed streams. Dataset/CLI
tests cover missing/duplicate/unknown options, complete integer parsing, invalid
labels and records, blank lines, long records, paths, and I/O errors.

First establish valid legacy characterization. For each defect or new contract,
run its test in RED, implement, then record GREEN. Retain commands and outcomes
in `docs/VALIDATION.md`; report measurements only when actually collected.

Integration tests train/save/reload/predict with actual upstream programs. Cover
the supplied regression files and an independently specified classification
fixture with expected held-out labels. Verify predictions and ordering, not only
process exit status or output existence. Do not claim generalization performance
from these correctness fixtures.

CI covers Linux/GCC, macOS/Clang, Windows/MSVC, Linux ASan/UBSan, project formatting,
dependency integrity, installation, and project-only coverage. Establish and enforce
a coverage floor from the measured baseline. Run the supported README workflow as
part of integration testing. All checks must pass on the PR's current commit before
merge; close the issue with evidence and remove the implementation worktree.

## Historical analyzer and future work

Keep `legacy/paired_sequence_analyzer.cpp` byte-for-byte unchanged and outside supported targets.
Restoration needs a paired-sequence data contract, orthogonal-array dependency
provenance, fixed-fold ownership, and a new public-API-only solver integration.
It is not made buildable by merely supplying upstream `svm.h`.

Python bindings, FASTA parsing, tuning/cross-validation, automatic alignment,
parallel execution, and model-quality research remain separately scoped features.
Their absence does not change the supported converter-to-official-solver workflow.
