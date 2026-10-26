"""Real solver, independent matrix oracle, and installed-artifact checks."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from itertools import combinations


def records(path: Path) -> list[tuple[str, float]]:
    return [(parts[0], float(parts[1])) for line in path.read_text(encoding="utf-8").splitlines()
            if (parts := line.split())]


def feature_map(sequence: str, kernel: str, k: int, decay: float) -> dict[object, float]:
    if kernel == "spectrum":
        features: dict[object, float] = {}
        for start in range(max(0, len(sequence) - k + 1)):
            token = sequence[start:start + k]
            features[token] = features.get(token, 0.0) + 1.0
        return features
    if kernel == "subsequence":
        features = {}
        for positions in combinations(range(len(sequence)), k):
            token = "".join(sequence[position] for position in positions)
            gaps = positions[-1] - positions[0] + 1 - k
            features[token] = features.get(token, 0.0) + decay ** gaps
        return features
    if kernel == "positional":
        return {(position, sequence[position]): 1.0 for position in range(len(sequence))}
    return {
        (position, order, sequence[position:position + order]):
            math.sqrt(2.0 * (k - order + 1) / (k * (k + 1)))
        for order in range(1, min(k, len(sequence)) + 1)
        for position in range(len(sequence) - order + 1)
    }


def oracle_rows(reference: list[tuple[str, float]], target: list[tuple[str, float]],
                kernel: str = "weighted-degree", k: int = 5,
                decay: float = 0.5) -> list[list[float]]:
    reference_features = [feature_map(sequence, kernel, k, decay) for sequence, _ in reference]
    result = []
    for sequence, _ in target:
        features = feature_map(sequence, kernel, k, decay)
        self_value = sum(value * value for value in features.values())
        row = []
        for other in reference_features:
            other_self = sum(value * value for value in other.values())
            raw = sum(value * other.get(key, 0.0) for key, value in features.items())
            row.append(raw / math.sqrt(self_value * other_self) if self_value and other_self else 0.0)
        result.append(row)
    return result


class LibsvmIntegrationTests(unittest.TestCase):
    options: argparse.Namespace

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="stringkernel-integration-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def run_program(self, *arguments: object) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [str(argument) for argument in arguments], capture_output=True, text=True,
            encoding="utf-8", timeout=45, check=False,
        )
        self.assertEqual(result.returncode, 0, f"{arguments!r}\n{result.stdout}\n{result.stderr}")
        return result

    def convert(self, reference: Path, target: Path | None, output: Path,
                converter: Path | None = None, kernel: str | None = None,
                k: int | None = None, decay: float | None = None) -> None:
        arguments = [converter or self.options.converter, "-train", reference]
        if target is not None:
            arguments.extend(["-test", target])
        if kernel is not None:
            arguments.extend(["--kernel", kernel])
        if k is not None:
            arguments.extend(["-k", str(k)])
        if decay is not None:
            arguments.extend(["--decay", str(decay)])
        result = self.run_program(*arguments)
        self.assertEqual(result.stderr, "")
        output.write_text(result.stdout, encoding="utf-8")

    def assert_matrix(self, path: Path, reference: list[tuple[str, float]],
                      target: list[tuple[str, float]], kernel: str = "weighted-degree",
                      k: int = 5, decay: float = 0.5) -> None:
        rows = path.read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(rows), len(target))
        for row_index, (line, expected) in enumerate(
                zip(rows, oracle_rows(reference, target, kernel, k, decay))):
            tokens = line.split()
            self.assertEqual(float(tokens[0]), target[row_index][1])
            self.assertEqual(len(tokens), len(reference) + 2)
            self.assertEqual(tokens[1], f"0:{row_index + 1}")
            for column, token in enumerate(tokens[2:]):
                index, value = token.split(":")
                self.assertEqual(int(index), column + 1)
                self.assertAlmostEqual(float(value), expected[column], delta=1e-12)

    def predictions(self, path: Path) -> list[float]:
        values = [float(line) for line in path.read_text(encoding="utf-8").splitlines()]
        self.assertTrue(all(math.isfinite(value) for value in values))
        return values

    def test_classification_persists_reloads_and_predicts_held_out_labels(self) -> None:
        source = self.options.source / "tests" / "fixtures"
        reference = source / "classification_train.txt"
        target = source / "classification_test.txt"
        training_matrix = self.directory / "classification.train"
        target_matrix = self.directory / "classification.test"
        model = self.directory / "classification.model"
        predictions = self.directory / "classification.predictions"
        self.convert(reference, None, training_matrix)
        self.convert(reference, target, target_matrix)
        self.assert_matrix(training_matrix, records(reference), records(reference))
        self.assert_matrix(target_matrix, records(reference), records(target))
        self.run_program(self.options.trainer, "-q", "-t", "4", training_matrix, model)
        self.assertIn("kernel_type precomputed\n", model.read_text(encoding="utf-8"))
        self.run_program(self.options.predictor, target_matrix, model, predictions)
        expected = [label for _, label in records(target)]
        self.assertEqual(self.predictions(predictions), expected)

    def test_new_kernels_match_oracles_and_predict_documented_examples(self) -> None:
        source = self.options.source / "examples" / "string-kernels"
        reference_path = source / "spectrum-train.tsv"
        target_path = source / "spectrum-test.tsv"
        reference = records(reference_path)
        target = records(target_path)
        expected = [label for _, label in target]
        for kernel, decay in [("spectrum", None), ("subsequence", 0.5),
                              ("subsequence", 0.25)]:
            with self.subTest(kernel=kernel, decay=decay):
                suffix = kernel + ("-quarter" if decay == 0.25 else "-default")
                training_matrix = self.directory / f"{suffix}.train"
                target_matrix = self.directory / f"{suffix}.test"
                model = self.directory / f"{suffix}.model"
                prediction_file = self.directory / f"{suffix}.predictions"
                self.convert(reference_path, None, training_matrix, kernel=kernel, k=2,
                             decay=decay)
                self.convert(reference_path, target_path, target_matrix, kernel=kernel, k=2,
                             decay=decay)
                parameter = 0.5 if decay is None else decay
                self.assert_matrix(training_matrix, reference, reference, kernel, 2, parameter)
                self.assert_matrix(target_matrix, reference, target, kernel, 2, parameter)
                self.run_program(self.options.trainer, "-q", "-t", "4", training_matrix, model)
                self.run_program(self.options.predictor, target_matrix, model, prediction_file)
                self.assertEqual(self.predictions(prediction_file), expected)

    def write_oracle_matrix(self, reference: list[tuple[str, float]],
                            target: list[tuple[str, float]], path: Path) -> None:
        lines = [f"{label!r} 0:{i + 1} " +
                 " ".join(f"{j + 1}:{value!r}" for j, value in enumerate(row))
                 for i, ((_, label), row) in enumerate(zip(target, oracle_rows(reference, target)))]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def test_supplied_regression_workflow_matches_independent_oracle_predictions(self) -> None:
        reference = self.options.source / "examples" / "regression" / "train_sequences.tsv"
        target = self.options.source / "examples" / "regression" / "test_sequences.tsv"
        generated_train = self.directory / "regression.train"
        generated_test = self.directory / "regression.test"
        oracle_train = self.directory / "oracle.train"
        oracle_test = self.directory / "oracle.test"
        self.convert(reference, None, generated_train)
        self.convert(reference, target, generated_test)
        self.assert_matrix(generated_train, records(reference), records(reference))
        self.assert_matrix(generated_test, records(reference), records(target))
        self.write_oracle_matrix(records(reference), records(reference), oracle_train)
        self.write_oracle_matrix(records(reference), records(target), oracle_test)
        outputs = []
        for prefix, training, testing in [("generated", generated_train, generated_test),
                                          ("oracle", oracle_train, oracle_test)]:
            model = self.directory / f"{prefix}.model"
            output = self.directory / f"{prefix}.predictions"
            self.run_program(self.options.trainer, "-q", "-s", "3", "-t", "4", training, model)
            self.assertIn("svm_type epsilon_svr\n", model.read_text(encoding="utf-8"))
            self.run_program(self.options.predictor, testing, model, output)
            outputs.append(self.predictions(output))
        self.assertEqual(len(outputs[0]), len(records(target)))
        self.assertEqual(len(outputs[1]), len(records(target)))
        for actual, expected in zip(*outputs):
            self.assertAlmostEqual(actual, expected, delta=1e-6)
        repeated = self.directory / "reloaded.predictions"
        self.run_program(self.options.predictor, generated_test,
                         self.directory / "generated.model", repeated)
        self.assertEqual(self.predictions(repeated), outputs[0])

    def test_installation_contains_public_artifacts_and_runs_the_documented_workflow(self) -> None:
        prefix = self.directory / "installed"
        self.run_program(self.options.cmake, "--install", self.options.build,
                         "--prefix", prefix, "--config", self.options.config)
        extension = self.options.converter.suffix
        converter = prefix / "bin" / f"StringKernelConverter{extension}"
        trainer = prefix / "bin" / f"svm-train{extension}"
        predictor = prefix / "bin" / f"svm-predict{extension}"
        for executable in [converter, trainer, predictor]:
            self.assertTrue(executable.is_file(), executable)
        for name in ["StringKernels.h", "Dataset.h", "LibsvmWriter.h"]:
            self.assertTrue((prefix / "include" / "StringKernel4LIBSVM" / name).is_file())
        self.assertTrue(list(prefix.rglob("libstringkernel.a")) or
                        list(prefix.rglob("stringkernel.lib")))
        copyright_file = prefix / "share" / "StringKernel4LIBSVM" / "COPYRIGHT"
        self.assertEqual(copyright_file.read_bytes(),
                         (self.options.source / "third_party" / "libsvm" / "COPYRIGHT").read_bytes())
        self.assertTrue((copyright_file.parent / "THIRD_PARTY_NOTICES.md").is_file())
        self.assertEqual(self.run_program(converter, "--version").stdout,
                         "StringKernel4LIBSVM 1.0.0 (LIBSVM 3.37)\n")
        train_matrix = self.directory / "installed.train"
        test_matrix = self.directory / "installed.test"
        model = self.directory / "installed.model"
        output = self.directory / "installed.predictions"
        example = self.options.source / "examples" / "regression"
        reference = example / "train_sequences.tsv"
        target = example / "test_sequences.tsv"
        self.convert(reference, None, train_matrix, converter)
        self.convert(reference, target, test_matrix, converter)
        self.run_program(trainer, "-q", "-s", "3", "-t", "4", train_matrix, model)
        self.run_program(predictor, test_matrix, model, output)
        self.assertEqual(len(self.predictions(output)), len(records(target)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for argument in ["converter", "trainer", "predictor", "source", "cmake", "build"]:
        parser.add_argument(f"--{argument}", type=Path, required=True)
    parser.add_argument("--config", default="Debug")
    options, remaining = parser.parse_known_args()
    for argument in ["converter", "trainer", "predictor", "source", "cmake", "build"]:
        setattr(options, argument, getattr(options, argument).resolve())
    LibsvmIntegrationTests.options = options
    unittest.main(argv=[__file__, *remaining])
