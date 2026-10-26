"""Process-level contracts for the supported converter CLI (standard library only)."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import subprocess
import tempfile
import unittest


class ConverterTests(unittest.TestCase):
    converter: Path

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="stringkernel-cli-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.train = self.directory / "training data.txt"
        self.train.write_text("AC 1\nAD -1\nZZ 2\n", encoding="utf-8")

    def run_converter(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(self.converter), *arguments], capture_output=True, text=True,
            encoding="utf-8", timeout=20, check=False,
        )

    def assert_failure(self, *arguments: str) -> None:
        result = self.run_converter(*arguments)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertEqual(result.stdout, "")
        self.assertTrue(result.stderr.strip())

    def test_valid_legacy_flags_and_analytical_kernel(self) -> None:
        result = self.run_converter("-train", str(self.train), "-k", "2", "-t", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        rows = result.stdout.splitlines()
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0].split()[1], "0:1")
        self.assertAlmostEqual(float(rows[0].split()[3].split(":")[1]), 0.4)
        self.assertEqual(float(rows[0].split()[4].split(":")[1]), 0)

    def test_help_and_pinned_versions(self) -> None:
        help_result = self.run_converter("--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("-train", help_result.stdout)
        self.assertIn("-test", help_result.stdout)
        version = self.run_converter("--version")
        self.assertEqual(version.returncode, 0, version.stderr)
        self.assertEqual(version.stdout, "StringKernel4LIBSVM 1.0.0 (LIBSVM 3.37)\n")

    def test_no_arguments_missing_file_and_missing_train_fail(self) -> None:
        self.assert_failure()
        self.assert_failure("-train", str(self.directory / "absent.txt"))
        self.assert_failure("-test", str(self.train))
        self.assert_failure("-k", "2")

    def test_missing_values_unknown_and_duplicate_options_fail(self) -> None:
        for arguments in [
            ["-train"], ["-train", str(self.train), "-test"],
            ["-train", str(self.train), "-k"], ["-train", str(self.train), "-t"],
            ["-train", str(self.train), "--wat"],
            ["-train", str(self.train), "-train", str(self.train)],
            ["--help", "-train", str(self.train)],
        ]:
            with self.subTest(arguments=arguments):
                self.assert_failure(*arguments)

    def test_invalid_integer_options_fail(self) -> None:
        for degree in ["0", "-1", "5.5", "x5", "5x", "", "2147483648"]:
            with self.subTest(degree=degree):
                self.assert_failure("-train", str(self.train), "-k", degree)
        for kernel in ["-1", "4", "1.0", "x1"]:
            with self.subTest(kernel=kernel):
                self.assert_failure("-train", str(self.train), "-t", kernel)

    def test_named_new_kernels_and_parameter_validation(self) -> None:
        self.train.write_text("ACG 1\nACT 1\nGGA -1\nGGT -1\n", encoding="utf-8")
        for kernel in ["spectrum", "subsequence"]:
            result = self.run_converter("-train", str(self.train), "--kernel", kernel, "-k", "2")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(result.stdout.splitlines()), 4)
        invalid = [
            ["--kernel", "spectrum", "-t", "2"],
            ["--kernel", "spectrum", "--decay", "0.5"],
            ["--kernel", "subsequence", "--decay", "nan"],
            ["--kernel", "subsequence", "--decay", "0"],
        ]
        for options in invalid:
            with self.subTest(options=options):
                self.assert_failure("-train", str(self.train), *options)

    def test_invalid_record_is_diagnosed_before_output(self) -> None:
        for record in ["BAD", "BAD nan", "BAD inf", "BAD 1e9999", "BAD 1 extra", "BAD 2x",
                       "BAD 0x10", "BAD +-1", "BAD ++1"]:
            with self.subTest(record=record):
                self.train.write_text(f"AC 1\n{record}\n", encoding="utf-8")
                result = self.run_converter("-train", str(self.train))
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(result.stdout, "")
                self.assertIn("2", result.stderr)
                self.assertIn(str(self.train), result.stderr)

    def test_empty_and_whitespace_only_datasets_fail(self) -> None:
        for contents in ["", " \t\r\n\n"]:
            self.train.write_text(contents, encoding="utf-8")
            self.assert_failure("-train", str(self.train))

    def test_target_errors_are_reported_before_output(self) -> None:
        target = self.directory / "target.txt"
        self.assert_failure("-train", str(self.train), "-test", str(target))
        target.write_text("AA 1\nBROKEN\n", encoding="utf-8")
        self.assert_failure("-train", str(self.train), "-test", str(target))

    def test_crlf_blanks_no_final_newline_and_case_sensitive_bytes(self) -> None:
        self.train.write_bytes(b"\r\nAC 1\r\n \t\r\nac 2")
        result = self.run_converter("-train", str(self.train))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 2)
        self.assertEqual(float(result.stdout.splitlines()[0].split()[-1].split(":")[1]), 0)

    def test_long_records_and_large_degree_are_bounded_by_sequence_length(self) -> None:
        self.train.write_text(
            "A" * 70000 + " 1\n" + "A" * 65534 + "B" * 4466 + " 2\n", encoding="utf-8"
        )
        result = self.run_converter("-train", str(self.train), "-k", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertAlmostEqual(float(result.stdout.splitlines()[0].split()[3].split(":")[1]),
                               65534 / 70000)
        self.train.write_text("AC 1\nAD -1\n", encoding="utf-8")
        large = self.run_converter("-train", str(self.train), "-k", "2147483647")
        self.assertEqual(large.returncode, 0, large.stderr)
        for row in large.stdout.splitlines():
            for token in row.split()[2:]:
                self.assertTrue(math.isfinite(float(token.split(":")[1])))

    def test_degree_five_is_default_and_unequal_wd_lengths_work(self) -> None:
        self.train.write_text("AAAA 1\nAAA 2\n", encoding="utf-8")
        default = self.run_converter("-train", str(self.train))
        explicit = self.run_converter("-train", str(self.train), "-k", "5", "-t", "1")
        self.assertEqual(default.returncode, 0, default.stderr)
        self.assertEqual(default.stdout, explicit.stdout)

    def test_dsk_requires_equal_lengths_in_reference_and_target(self) -> None:
        self.train.write_text("AA 1\nA 2\n", encoding="utf-8")
        self.assert_failure("-train", str(self.train), "-t", "0")
        self.train.write_text("AA 1\nAB 2\n", encoding="utf-8")
        target = self.directory / "target.txt"
        target.write_text("A 1\n", encoding="utf-8")
        self.assert_failure("-train", str(self.train), "-test", str(target), "-t", "0")

    def test_output_preserves_regression_label_precision(self) -> None:
        for value in ["0.12345678901234566", "4.9406564584124654e-324",
                      "1.7976931348623157e308", "+1.25"]:
            with self.subTest(value=value):
                self.train.write_text(f"AC {value}\n", encoding="utf-8")
                result = self.run_converter("-train", str(self.train))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(float(result.stdout.split()[0]), float(value))

    @unittest.skipUnless(Path("/dev/full").exists(), "requires Linux /dev/full")
    def test_redirected_output_failure_has_nonzero_status(self) -> None:
        for arguments in [["-train", str(self.train)], ["--help"], ["--version"]]:
            with self.subTest(arguments=arguments):
                with Path("/dev/full").open("w", encoding="utf-8") as output:
                    result = subprocess.run(
                        [str(self.converter), *arguments], stdout=output,
                        stderr=subprocess.PIPE, text=True, timeout=20, check=False,
                    )
                self.assertNotEqual(result.returncode, 0)
                self.assertTrue(result.stderr.strip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--converter", type=Path, required=True)
    options, remaining = parser.parse_known_args()
    ConverterTests.converter = options.converter.resolve()
    unittest.main(argv=[__file__, *remaining])
