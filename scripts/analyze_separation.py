#!/usr/bin/env python3
"""Render deterministic, fixture-specific string-kernel separation examples."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from html import escape
from itertools import combinations
import math
from pathlib import Path
import subprocess
import tempfile


@dataclass(frozen=True)
class Record:
    sequence: str
    label: int


@dataclass(frozen=True)
class Example:
    name: str
    filename: str
    kernel: str
    degree: int
    decay: float
    lesson: str


EXAMPLES = (
    Example("aligned motif", "aligned-motif.tsv", "weighted-degree", 3, 0.5,
            "The class cue stays at the same positions, so aligned local matches help."),
    Example("moved motif", "moved-motif.tsv", "spectrum", 3, 0.5,
            "The same exact 3-mer occurs at different positions, so position-free counts help."),
    Example("gapped motif", "gapped-motif.tsv", "subsequence", 3, 0.5,
            "The ordered class cue survives inserted bytes, so a gap-tolerant subsequence helps."),
)


def read_records(path: Path) -> tuple[list[Record], list[Record]]:
    training: list[Record] = []
    testing: list[Record] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split()
        if len(fields) != 3 or fields[0] not in {"train", "test"}:
            raise ValueError(f"{path}:{line_number}: expected 'train|test SEQUENCE LABEL'")
        sequence, label_token = fields[1:]
        try:
            label = int(label_token)
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: label must be -1 or 1") from error
        if not sequence or label not in {-1, 1}:
            raise ValueError(f"{path}:{line_number}: expected a sequence and label -1 or 1")
        record = Record(sequence, label)
        (training if fields[0] == "train" else testing).append(record)
    if len(training) < 2 or len(testing) != 2 or {r.label for r in training} != {-1, 1}:
        raise ValueError(f"{path}: expected both training classes and two held-out records")
    if [r.label for r in testing] != [1, -1]:
        raise ValueError(f"{path}: held-out labels must be ordered positive then negative")
    return training, testing


def feature_map(sequence: str, example: Example) -> dict[object, float]:
    k = example.degree
    if example.kernel == "weighted-degree":
        features: dict[object, float] = {}
        for order in range(1, min(k, len(sequence)) + 1):
            beta = 2.0 * (k - order + 1) / (k * (k + 1))
            scale = math.sqrt(beta)
            for start in range(len(sequence) - order + 1):
                features[(start, order, sequence[start:start + order])] = scale
        return features
    if example.kernel == "spectrum":
        features = {}
        for start in range(max(0, len(sequence) - k + 1)):
            token = sequence[start:start + k]
            features[token] = features.get(token, 0.0) + 1.0
        return features
    features = {}
    for positions in combinations(range(len(sequence)), k):
        token = "".join(sequence[position] for position in positions)
        gap_count = positions[-1] - positions[0] + 1 - k
        features[token] = features.get(token, 0.0) + example.decay ** gap_count
    return features


def normalized(left: dict[object, float], right: dict[object, float]) -> float:
    left_self = sum(value * value for value in left.values())
    right_self = sum(value * value for value in right.values())
    if not left_self or not right_self:
        return 0.0
    raw = sum(value * right.get(key, 0.0) for key, value in left.items())
    return raw / math.sqrt(left_self * right_self)


def expected_matrix(reference: list[Record], targets: list[Record],
                    example: Example) -> list[list[float]]:
    reference_maps = [feature_map(record.sequence, example) for record in reference]
    return [[normalized(feature_map(target.sequence, example), ref_map)
             for ref_map in reference_maps] for target in targets]


def parse_matrix(text: str, expected_records: list[Record],
                 expected_columns: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for row_index, line in enumerate(text.splitlines(), 1):
        fields = line.split()
        if len(fields) != expected_columns + 2:
            raise ValueError("converter emitted an unexpected matrix width")
        if float(fields[0]) != expected_records[row_index - 1].label:
            raise ValueError("converter changed sequence-label order")
        if fields[1] != f"0:{row_index}":
            raise ValueError("converter emitted an unexpected row serial")
        values: list[float] = []
        for column, field in enumerate(fields[2:], 1):
            index, value = field.split(":", 1)
            if index != str(column):
                raise ValueError("converter changed reference column order")
            values.append(float(value))
        rows.append(values)
    if len(rows) != len(expected_records):
        raise ValueError("converter emitted an unexpected matrix height")
    return rows


def run(arguments: list[object]) -> str:
    result = subprocess.run([str(value) for value in arguments],
                            capture_output=True, text=True, encoding="utf-8", timeout=45)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {arguments!r}\n{result.stderr}")
    return result.stdout


def convert(converter: Path, train_path: Path, test_path: Path | None,
            example: Example) -> str:
    command: list[object] = [converter, "-train", train_path,
                             "--kernel", example.kernel, "-k", example.degree]
    if test_path is not None:
        command.extend(["-test", test_path])
    if example.kernel == "subsequence":
        command.extend(["--decay", example.decay])
    return run(command)


def class_similarity(training: list[Record], matrix: list[list[float]]) -> tuple[float, float]:
    within: list[float] = []
    between: list[float] = []
    for left in range(len(training)):
        for right in range(left + 1, len(training)):
            (within if training[left].label == training[right].label else between).append(
                matrix[left][right])
    return sum(within) / len(within), sum(between) / len(between)


def fmt(value: float) -> str:
    return f"{value:.3f}"


def color(value: float) -> str:
    value = min(1.0, max(0.0, value))
    start = (249, 246, 255)
    end = (72, 35, 116)
    rgb = tuple(round(a + (b - a) * value) for a, b in zip(start, end))
    return "#" + "".join(f"{component:02x}" for component in rgb)


def svg_panel(example: Example, records: list[Record], matrix: list[list[float]], x: int) -> str:
    cell = 20
    left = x + 74
    top = 61
    parts = [f'<g><text x="{x+4}" y="22" class="title">{example.name}</text>',
             f'<text x="{x+4}" y="39" class="subtitle">{example.kernel} · k={example.degree}'
             f'{f" · λ={example.decay}" if example.kernel == "subsequence" else ""}</text>']
    for index, record in enumerate(records):
        y = top + index * cell
        parts.append(f'<text x="{x+4}" y="{y+14}" class="label">'
                     f'{escape(record.sequence)} {"+" if record.label == 1 else "−"}</text>')
        for column, reference in enumerate(records[:len(matrix[0])]):
            value = matrix[index][column]
            cell_x = left + column * cell
            parts.append(f'<rect x="{cell_x}" y="{y}" width="{cell}" height="{cell}" '
                         f'fill="{color(value)}" stroke="#ffffff"><title>'
                         f'{escape(record.sequence)} vs {escape(reference.sequence)}: '
                         f'{value:.6f}</title></rect>')
    for index, record in enumerate(records[:len(matrix[0])]):
        angle_x = left + index * cell + 3
        parts.append(f'<text transform="translate({angle_x},{top+len(records)*cell+6}) rotate(55)" '
                     f'class="axis">{escape(record.sequence)}</text>')
    parts.append("</g>")
    return "".join(parts)


def render_svg(panels: list[tuple[Example, list[Record], list[list[float]]]]) -> str:
    chunks = ["<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"780\" height=\"270\" "
              "viewBox=\"0 0 780 270\"><style>text{font-family:system-ui,sans-serif}" 
              ".title{font-size:13px;font-weight:700}.subtitle{font-size:10px;fill:#555}" 
              ".label,.axis{font-size:9px;fill:#333}</style>"]
    for index, panel in enumerate(panels):
        chunks.append(svg_panel(*panel, x=index * 260))
    chunks.append("<text x=\"8\" y=\"246\" class=\"subtitle\">Darker = more similar; "
                  "rows are examples, columns are training references.</text></svg>\n")
    return "".join(chunks)


def analyze(converter: Path, trainer: Path, predictor: Path,
            source: Path) -> tuple[str, str]:
    panels: list[tuple[Example, list[Record], list[list[float]]]] = []
    lines = ["# A small visual tour of string-kernel separation", "",
             "These are **fixture-specific** geometry demonstrations, not representative "
             "datasets or estimates of real-world accuracy.", "",
             "Each heatmap combines training rows and held-out rows against the original "
             "training references. Darker cells mean higher normalized similarity.", ""]
    with tempfile.TemporaryDirectory(prefix="stringkernel-analysis-data-") as temporary:
        scratch = Path(temporary)
        for example in EXAMPLES:
            training, testing = read_records(source / "examples" / "analysis" / example.filename)
            train_path = scratch / f"{example.filename}.train"
            test_path = scratch / f"{example.filename}.test"
            train_path.write_text("".join(f"{row.sequence} {row.label}\n" for row in training),
                                  encoding="utf-8")
            test_path.write_text("".join(f"{row.sequence} {row.label}\n" for row in testing),
                                 encoding="utf-8")
            train_matrix_path = scratch / f"{example.filename}.matrix"
            test_matrix_path = scratch / f"{example.filename}.target"
            train_text = convert(converter, train_path, None, example)
            test_text = convert(converter, train_path, test_path, example)
            train_matrix_path.write_text(train_text, encoding="utf-8")
            test_matrix_path.write_text(test_text, encoding="utf-8")
            train_values = parse_matrix(train_text, training, len(training))
            test_values = parse_matrix(test_text, testing, len(training))
            oracle_train = expected_matrix(training, training, example)
            oracle_test = expected_matrix(training, testing, example)
            for actual_rows, oracle_rows in [(train_values, oracle_train),
                                             (test_values, oracle_test)]:
                for actual, expected in zip(actual_rows, oracle_rows):
                    if any(not math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12)
                           for a, b in zip(actual, expected)):
                        raise ValueError(f"{example.name}: converter disagrees with "
                                         "feature-map oracle")
            model_path = scratch / f"{example.filename}.model"
            prediction_path = scratch / f"{example.filename}.predictions"
            run([trainer, "-q", "-t", "4", train_matrix_path, model_path])
            # LIBSVM prediction consumes the test-to-training precomputed matrix.
            run([predictor, test_matrix_path, model_path, prediction_path])
            predictions = [float(line) for line in
                           prediction_path.read_text(encoding="utf-8").splitlines()]
            expected_labels = [record.label for record in testing]
            if predictions != [float(label) for label in expected_labels]:
                raise ValueError(f"{example.name}: expected held-out labels "
                                 f"{expected_labels}, got {predictions}")
            combined_records = training + testing
            combined_matrix = train_values + test_values
            within, between = class_similarity(training, train_values)
            lines.extend([f"## {example.name}", "",
                          f"**Selected kernel:** `{example.kernel}`" +
                          (f" (`k={example.degree}`, `λ={example.decay}`)"
                           if example.kernel == "subsequence"
                           else f" (`k={example.degree}`)"), "",
                          example.lesson, "",
                          f"Training-pair mean similarity: within class **{fmt(within)}**; "
                          f"between classes **{fmt(between)}**. These values summarize only "
                          "the four constructed training strings.", "",
                          f"LIBSVM held-out predictions: "
                          f"`{', '.join(str(int(value)) for value in predictions)}` "
                          f"(expected `{', '.join(str(value) for value in expected_labels)}`).",
                          ""])
            panels.append((example, combined_records, combined_matrix))
    lines.extend(["## Reading the picture", "",
                  "The aligned example rewards same-position local matches; the moved example "
                  "rewards exact k-mers wherever they occur; the gapped example rewards "
                  "ordered patterns despite inserted bytes. Each panel uses the kernel "
                  "chosen for that construction, so this is a demonstration of fit rather "
                  "than a head-to-head scorecard.", "",
                  "All records are synthetic and intentionally tiny. The labels and test "
                  "examples were designed together with the patterns. The predictions "
                  "verify the converter-to-LIBSVM workflow; they do not measure "
                  "generalization, biological usefulness, or comparative runtime.", ""])
    return "\n".join(lines), render_svg(panels)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--converter", type=Path, required=True)
    parser.add_argument("--trainer", type=Path, required=True)
    parser.add_argument("--predictor", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    try:
        relative_output = output.relative_to(source)
    except ValueError:
        relative_output = None
    if relative_output is not None and (
            not relative_output.parts or relative_output.parts[0] != "build"):
        parser.error("output directory must be outside the source checkout or "
                     "inside its build directory")
    try:
        markdown, graphic = analyze(args.converter, args.trainer, args.predictor, source)
        output.mkdir(parents=True, exist_ok=True)
        with (output / "analysis.md").open("w", encoding="utf-8", newline="\n") as report:
            report.write(markdown)
        with (output / "similarities.svg").open("w", encoding="utf-8", newline="\n") as image:
            image.write(graphic)
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
        parser.error(str(error))
    print(f"Wrote {output / 'analysis.md'} and {output / 'similarities.svg'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
