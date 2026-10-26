"""End-to-end checks for the reproducible kernel-separation report."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from xml.etree import ElementTree


class SeparationAnalysisTests(unittest.TestCase):
    options: argparse.Namespace

    def test_report_is_reproducible_and_models_recover_fixture_labels(self) -> None:
        with tempfile.TemporaryDirectory(prefix="stringkernel-analysis-") as temporary:
            output = Path(temporary) / "report"
            command = [
                sys.executable, str(self.options.source / "scripts" / "analyze_separation.py"),
                "--converter", str(self.options.converter),
                "--trainer", str(self.options.trainer),
                "--predictor", str(self.options.predictor),
                "--output-dir", str(output),
            ]
            first = subprocess.run(command, capture_output=True, text=True, timeout=90)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            report = (output / "analysis.md").read_bytes()
            graphic = (output / "similarities.svg").read_bytes()
            self.assertIn(b"fixture-specific", report)
            self.assertIn(b"aligned motif", report)
            self.assertIn(b"moved motif", report)
            self.assertIn(b"gapped motif", report)
            self.assertEqual(report.count(b"LIBSVM held-out predictions: `1, -1`"), 3)
            self.assertIn(b"within class", report)
            self.assertIn(b"between classes", report)
            svg = ElementTree.fromstring(graphic)
            self.assertEqual(svg.tag, "{http://www.w3.org/2000/svg}svg")
            self.assertIn(b"spectrum", graphic)
            self.assertIn(b"subsequence", graphic)
            published_figure = self.options.source / "docs" / "analysis" / "similarities.svg"
            self.assertEqual(graphic, published_figure.read_bytes())
            first_paths = sorted(path.relative_to(output).as_posix()
                                 for path in output.rglob("*") if path.is_file())

            second = subprocess.run(command, capture_output=True, text=True, timeout=90)
            self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
            self.assertEqual((output / "analysis.md").read_bytes(), report)
            self.assertEqual((output / "similarities.svg").read_bytes(), graphic)
            self.assertEqual(sorted(path.relative_to(output).as_posix()
                                    for path in output.rglob("*") if path.is_file()), first_paths)
            self.assertEqual(first_paths, ["analysis.md", "similarities.svg"])

            readme = self.options.source / "README.md"
            original_readme = readme.read_bytes()
            command[command.index("--output-dir") + 1] = str(readme)
            rejected = subprocess.run(command, capture_output=True, text=True, timeout=90)
            self.assertEqual(rejected.returncode, 2, rejected.stdout + rejected.stderr)
            self.assertIn("outside the source checkout", rejected.stderr)
            self.assertEqual(readme.read_bytes(), original_readme)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--converter", type=Path, required=True)
    parser.add_argument("--trainer", type=Path, required=True)
    parser.add_argument("--predictor", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    options, remaining = parser.parse_known_args()
    SeparationAnalysisTests.options = options
    unittest.main(argv=[__file__, *remaining])


if __name__ == "__main__":
    main()
