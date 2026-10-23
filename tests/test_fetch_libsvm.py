"""Behavioral contracts for the standard-library LIBSVM bootstrap helper."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import copy
import hashlib
from http.client import IncompleteRead
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import URLError
from urllib.request import Request


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "fetch_libsvm.py"
COMMIT = "6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb"
REQUIRED_FILES = (
    "svm.cpp", "svm.h", "svm-train.c", "svm-predict.c", "README", "COPYRIGHT",
)
CONTENTS = {
    name: f"Fixture LIBSVM content for {name}\n".encode("utf-8")
    for name in REQUIRED_FILES
}


def fixture_manifest() -> dict:
    return {
        "name": "LIBSVM",
        "version": "3.37",
        "tag": "v337",
        "commit": COMMIT,
        "repository": "https://github.com/cjlin1/libsvm",
        "files": {
            name: {
                "url": f"https://raw.githubusercontent.com/cjlin1/libsvm/{COMMIT}/{name}",
                "sha256": hashlib.sha256(data).hexdigest(),
            }
            for name, data in CONTENTS.items()
        },
    }


class FetchLibsvmTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        spec = importlib.util.spec_from_file_location("fetch_libsvm_under_test", SCRIPT)
        if spec is None or spec.loader is None:
            raise AssertionError(f"cannot import bootstrap script: {SCRIPT}")
        cls.fetcher = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = cls.fetcher
        cls.addClassCleanup(sys.modules.pop, spec.name, None)
        spec.loader.exec_module(cls.fetcher)

    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="stringkernel-fetch-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.source = self.directory / "sources with spaces" / "libsvm"
        self.manifest = self.directory / "libsvm.json"
        self.metadata = fixture_manifest()
        self.write_manifest(self.metadata)

    def write_manifest(self, value: object) -> None:
        self.manifest.write_text(json.dumps(value), encoding="utf-8")

    def seed_cache(self, names=REQUIRED_FILES) -> None:
        self.source.mkdir(parents=True, exist_ok=True)
        for name in names:
            (self.source / name).write_bytes(CONTENTS[name])

    def assert_complete_cache(self) -> None:
        for name, expected in CONTENTS.items():
            with self.subTest(file=name):
                actual = (self.source / name).read_bytes()
                self.assertEqual(actual, expected)
                self.assertEqual(
                    hashlib.sha256(actual).hexdigest(),
                    self.metadata["files"][name]["sha256"],
                )

    def assert_unpublished(self, names=REQUIRED_FILES) -> None:
        for name in names:
            with self.subTest(file=name):
                self.assertFalse((self.source / name).exists(), f"published {name} on failure")

    def http_response(self, request, *args, **kwargs) -> io.BytesIO:
        url = request.full_url if isinstance(request, Request) else request
        responses = {
            self.metadata["files"][name]["url"]: data
            for name, data in CONTENTS.items()
        }
        self.assertIn(url, responses, "fetch must use the pinned official URL")
        return io.BytesIO(responses[url])

    @staticmethod
    def requested_urls(mock) -> list[str]:
        return [
            call.args[0].full_url if isinstance(call.args[0], Request) else call.args[0]
            for call in mock.call_args_list
        ]

    def assert_integrity_diagnostic(self, error: Exception, name: str) -> None:
        message = str(error)
        self.assertIn(name, message)
        self.assertRegex(message.lower(), r"sha-?256|digest|hash")

    def test_missing_cache_downloads_all_pinned_bytes(self) -> None:
        with patch("urllib.request.urlopen", side_effect=self.http_response) as download:
            self.fetcher.ensure_libsvm(self.source, self.manifest)
        self.assert_complete_cache()
        expected = [entry["url"] for entry in self.metadata["files"].values()]
        self.assertCountEqual(self.requested_urls(download), expected)

    def test_complete_valid_cache_works_offline(self) -> None:
        self.seed_cache()
        with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as download:
            self.fetcher.ensure_libsvm(self.source, self.manifest)
        download.assert_not_called()
        self.assert_complete_cache()

    def test_partial_valid_cache_fetches_only_missing_files_and_then_works_offline(self) -> None:
        existing = ("svm.cpp", "README")
        self.seed_cache(existing)
        with patch("urllib.request.urlopen", side_effect=self.http_response) as download:
            self.fetcher.ensure_libsvm(self.source, self.manifest)
        expected = [
            self.metadata["files"][name]["url"]
            for name in REQUIRED_FILES if name not in existing
        ]
        self.assertCountEqual(self.requested_urls(download), expected)
        self.assert_complete_cache()
        with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as offline:
            self.fetcher.ensure_libsvm(self.source, self.manifest)
        offline.assert_not_called()

    def test_existing_corruption_is_untouched_and_rejected_before_network(self) -> None:
        self.seed_cache(("svm.cpp", "COPYRIGHT"))
        corrupt = self.source / "COPYRIGHT"
        corrupt.write_bytes(b"locally modified upstream notice\n")
        before = {path.name: path.read_bytes() for path in self.source.iterdir()}
        with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as download:
            with self.assertRaises(self.fetcher.BootstrapError) as caught:
                self.fetcher.ensure_libsvm(self.source, self.manifest)
        download.assert_not_called()
        self.assert_integrity_diagnostic(caught.exception, "COPYRIGHT")
        self.assertRegex(str(caught.exception).lower(), r"remove|delete")
        self.assertEqual({path.name: path.read_bytes() for path in self.source.iterdir()}, before)

    def test_bad_download_hash_publishes_none_of_the_missing_files(self) -> None:
        self.seed_cache(("svm.h",))

        def response(request, *args, **kwargs):
            url = request.full_url if isinstance(request, Request) else request
            if url == self.metadata["files"]["COPYRIGHT"]["url"]:
                return io.BytesIO(b"wrong upstream content\n")
            return self.http_response(request, *args, **kwargs)

        with patch("urllib.request.urlopen", side_effect=response):
            with self.assertRaises(self.fetcher.BootstrapError) as caught:
                self.fetcher.ensure_libsvm(self.source, self.manifest)
        self.assert_integrity_diagnostic(caught.exception, "COPYRIGHT")
        self.assertEqual((self.source / "svm.h").read_bytes(), CONTENTS["svm.h"])
        self.assert_unpublished(name for name in REQUIRED_FILES if name != "svm.h")
        with patch("urllib.request.urlopen", side_effect=self.http_response):
            self.fetcher.ensure_libsvm(self.source, self.manifest)
        self.assert_complete_cache()

    def test_network_failure_is_actionable_and_does_not_publish_partial_files(self) -> None:
        def response(request, *args, **kwargs):
            url = request.full_url if isinstance(request, Request) else request
            if url == self.metadata["files"]["COPYRIGHT"]["url"]:
                raise URLError("fixture network unavailable")
            return self.http_response(request, *args, **kwargs)

        with patch("urllib.request.urlopen", side_effect=response):
            with self.assertRaises(self.fetcher.BootstrapError) as caught:
                self.fetcher.ensure_libsvm(self.source, self.manifest)
        message = str(caught.exception)
        self.assertIn("COPYRIGHT", message)
        self.assertIn("fixture network unavailable", message)
        self.assertRegex(message.lower(), r"retry|rerun|fetch_libsvm\.py")
        self.assert_unpublished()

    def test_interrupted_response_is_actionable_and_does_not_publish_partial_files(self) -> None:
        class InterruptedResponse(io.BytesIO):
            def read(self, *args, **kwargs):
                raise IncompleteRead(b"partial", 17)

        def response(request, *args, **kwargs):
            url = request.full_url if isinstance(request, Request) else request
            if url == self.metadata["files"]["COPYRIGHT"]["url"]:
                return InterruptedResponse()
            return self.http_response(request, *args, **kwargs)

        with patch("urllib.request.urlopen", side_effect=response):
            with self.assertRaises(self.fetcher.BootstrapError) as caught:
                self.fetcher.ensure_libsvm(self.source, self.manifest)
        message = str(caught.exception)
        self.assertIn("COPYRIGHT", message)
        self.assertRegex(message.lower(), r"incomplete|read|download")
        self.assertRegex(message.lower(), r"retry|rerun|fetch_libsvm\.py")
        self.assert_unpublished()

    def test_invalid_source_directory_is_diagnosed_without_network(self) -> None:
        self.source.parent.mkdir(parents=True)
        self.source.write_bytes(b"a file is not a source directory\n")
        with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as download:
            with self.assertRaises(self.fetcher.BootstrapError) as caught:
                self.fetcher.ensure_libsvm(self.source, self.manifest)
        download.assert_not_called()
        self.assertIn(str(self.source), str(caught.exception))
        self.assertIn("director", str(caught.exception).lower())

    def test_directory_in_place_of_a_required_file_is_rejected_without_network(self) -> None:
        self.source.mkdir(parents=True)
        (self.source / "svm.cpp").mkdir()
        with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as download:
            with self.assertRaises(self.fetcher.BootstrapError) as caught:
                self.fetcher.ensure_libsvm(self.source, self.manifest)
        download.assert_not_called()
        self.assertIn("svm.cpp", str(caught.exception))

    def assert_manifest_rejected(self, value: object) -> None:
        self.write_manifest(value)
        with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as download:
            with self.assertRaises(self.fetcher.BootstrapError) as caught:
                self.fetcher.ensure_libsvm(self.source, self.manifest)
        download.assert_not_called()
        self.assertIn("manifest", str(caught.exception).lower())
        self.assert_unpublished()

    def test_manifest_requires_the_exact_six_filenames(self) -> None:
        for replacement in (None, "unexpected.c", "../svm.cpp", "..\\svm.cpp", "/svm.cpp"):
            with self.subTest(replacement=replacement):
                value = copy.deepcopy(self.metadata)
                entry = value["files"].pop("svm.cpp")
                if replacement is not None:
                    value["files"][replacement] = entry
                self.assert_manifest_rejected(value)
        value = copy.deepcopy(self.metadata)
        value["files"]["unexpected.c"] = value["files"]["svm.cpp"]
        self.assert_manifest_rejected(value)

    def test_manifest_requires_a_lowercase_full_commit_id(self) -> None:
        for commit in ("master", COMMIT[:-1], COMMIT.upper(), "g" * 40, None, 337):
            with self.subTest(commit=commit):
                value = copy.deepcopy(self.metadata)
                value["commit"] = commit
                self.assert_manifest_rejected(value)

    def test_manifest_rejects_urls_that_are_not_the_exact_official_pin(self) -> None:
        for url in (
            f"https://example.com/cjlin1/libsvm/{COMMIT}/svm.cpp",
            "https://raw.githubusercontent.com/cjlin1/libsvm/master/svm.cpp",
            f"http://raw.githubusercontent.com/cjlin1/libsvm/{COMMIT}/svm.cpp",
            f"https://raw.githubusercontent.com/cjlin1/libsvm/{COMMIT}/svm.h",
            f"https://raw.githubusercontent.com/cjlin1/libsvm/{COMMIT}/svm.cpp?redirect=1",
            None,
        ):
            with self.subTest(url=url):
                value = copy.deepcopy(self.metadata)
                value["files"]["svm.cpp"]["url"] = url
                self.assert_manifest_rejected(value)

    def test_manifest_requires_sha256_digests(self) -> None:
        for digest in ("", "0" * 63, "g" * 64, None, 337):
            with self.subTest(digest=digest):
                value = copy.deepcopy(self.metadata)
                value["files"]["svm.cpp"]["sha256"] = digest
                self.assert_manifest_rejected(value)

    def test_manifest_requires_nonempty_version_and_tag(self) -> None:
        for key in ("version", "tag"):
            for invalid in ("", " \t", None, 337):
                with self.subTest(field=key, value=invalid):
                    value = copy.deepcopy(self.metadata)
                    value[key] = invalid
                    self.assert_manifest_rejected(value)

    def test_manifest_rejects_invalid_structure(self) -> None:
        invalid_values = [[], None, {}, {"files": []}]
        for files in (None, [], "svm.cpp"):
            value = copy.deepcopy(self.metadata)
            value["files"] = files
            invalid_values.append(value)
        for entry in (None, [], "not a file entry", {}):
            value = copy.deepcopy(self.metadata)
            value["files"]["svm.cpp"] = entry
            invalid_values.append(value)
        for value in invalid_values:
            with self.subTest(manifest=value):
                self.assert_manifest_rejected(value)

    def test_missing_or_malformed_manifest_is_actionable(self) -> None:
        for contents in (None, "{not JSON", "[1, 2,"):
            with self.subTest(contents=contents):
                if contents is None:
                    self.manifest.unlink()
                else:
                    self.manifest.write_text(contents, encoding="utf-8")
                with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as download:
                    with self.assertRaises(self.fetcher.BootstrapError) as caught:
                        self.fetcher.ensure_libsvm(self.source, self.manifest)
                download.assert_not_called()
                self.assertIn(str(self.manifest), str(caught.exception))
                self.assert_unpublished()

    def test_main_reports_invalid_source_directory_on_stderr(self) -> None:
        self.source.parent.mkdir(parents=True)
        self.source.write_text("not a directory", encoding="utf-8")
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network")) as download:
            with redirect_stdout(stdout), redirect_stderr(stderr):
                result = self.fetcher.main(["--source-dir", str(self.source)])
        self.assertIsInstance(result, int)
        self.assertNotEqual(result, 0)
        self.assertIn(str(self.source), stderr.getvalue())
        self.assertIn("director", stderr.getvalue().lower())
        self.assertEqual(stdout.getvalue(), "")
        download.assert_not_called()

    def test_main_reports_network_failure_on_stderr(self) -> None:
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch("urllib.request.urlopen", side_effect=URLError("fixture network unavailable")):
            with redirect_stdout(stdout), redirect_stderr(stderr):
                result = self.fetcher.main(["--source-dir", str(self.source)])
        self.assertIsInstance(result, int)
        self.assertNotEqual(result, 0)
        self.assertIn("fixture network unavailable", stderr.getvalue())
        self.assertRegex(stderr.getvalue().lower(), r"retry|rerun|fetch_libsvm\.py")
        self.assertEqual(stdout.getvalue(), "")
        self.assert_unpublished()


class FetchLibsvmCliTests(unittest.TestCase):
    def test_help_is_a_real_cli_command(self) -> None:
        with tempfile.TemporaryDirectory(prefix="stringkernel-fetch-cli-") as directory:
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--help"],
                cwd=directory, capture_output=True, text=True, encoding="utf-8",
                timeout=20, check=False,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--source-dir", result.stdout)
        self.assertEqual(result.stderr, "")


if __name__ == "__main__":
    unittest.main()
