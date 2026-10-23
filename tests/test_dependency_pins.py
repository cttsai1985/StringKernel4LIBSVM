"""Pin upstream content independently of its editable local provenance manifest."""

import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent.parent
COMMIT = "6b907139084abf2da4d6d3cb10dc3b7eaffa2fbb"
EXPECTED_HASHES = {
    "svm.cpp": "bfe02a7c001151ab703c206ff6a9136c75eb704b4427536e3f65d3ac90bd0eb6",
    "svm.h": "597aff50cfb9a04922736472bd125862439c87408a83cb57fb04980744dffd51",
    "svm-train.c": "004f0ba35101a8564fd535c13abacfc3e166f8ebddda73729f4be93526cc15bc",
    "svm-predict.c": "01bacc60534bfb07b22ad3283e084a31e795d2835981e743eefa42d048981194",
    "README": "6322d5062bedcb8148d5320588909325568fe72247df2d23c8c6be8b6d8c0743",
    "COPYRIGHT": "958fd47be5e5b94dabeb68b077fbe2adbbd49c015fe2e453fe11afeb4a73c7e5",
}
CATCH2_SHA256 = "d16bd1b8b2364918bd472cb216394acb93bde6784890df0d91bccfb3f28d777d"


class DependencyPins(unittest.TestCase):
    def test_libsvm_files_match_independent_pins(self):
        for name, expected_hash in EXPECTED_HASHES.items():
            with self.subTest(file=name):
                path = ROOT / "third_party" / "libsvm" / name
                self.assertTrue(path.is_file(), f"missing bootstrapped source: {name}")
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), expected_hash)

    def test_manifest_agrees_with_independent_pins(self):
        path = ROOT / "cmake" / "libsvm.json"
        self.assertTrue(path.is_file(), "missing upstream provenance")
        manifest = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["version"], "3.37")
        self.assertEqual(manifest["tag"], "v337")
        self.assertEqual(manifest["commit"], COMMIT)
        self.assertEqual(set(manifest["files"]), set(EXPECTED_HASHES))
        for name, expected_hash in EXPECTED_HASHES.items():
            with self.subTest(file=name):
                entry = manifest["files"][name]
                self.assertEqual(entry["sha256"], expected_hash)
                self.assertEqual(
                    entry["url"],
                    f"https://raw.githubusercontent.com/cjlin1/libsvm/{COMMIT}/{name}",
                )

    def check_cmake_cache(self, cache, expected_success):
        with tempfile.TemporaryDirectory(prefix="stringkernel-cmake-pins-") as directory:
            source = Path(directory)
            (source / "cmake").mkdir()
            shutil.copyfile(ROOT / "cmake" / "libsvm.json", source / "cmake" / "libsvm.json")
            if cache:
                shutil.copytree(ROOT / "third_party" / "libsvm", source / "third_party" / "libsvm")
            if cache == "corrupt":
                (source / "third_party" / "libsvm" / "README").write_bytes(b"modified")
            script = source / "verify.cmake"
            script.write_text(
                f'set(PROJECT_SOURCE_DIR "{source.as_posix()}")\n'
                f'include("{(ROOT / "cmake" / "Libsvm.cmake").as_posix()}")\n',
                encoding="utf-8",
            )
            result = subprocess.run(["cmake", "-P", str(script)], capture_output=True,
                                    text=True, timeout=30, check=False)
            output = result.stdout + result.stderr
            self.assertEqual(result.returncode == 0, expected_success, output)
            return output

    def test_cmake_rejects_missing_sources_with_bootstrap_command(self):
        self.assertIn("scripts/fetch_libsvm.py", self.check_cmake_cache(None, False))

    def test_cmake_rejects_modified_cached_source(self):
        output = self.check_cmake_cache("corrupt", False)
        self.assertIn("README", output)
        self.assertIn("SHA-256", output)
        self.assertIn("scripts/fetch_libsvm.py", output)

    def test_cmake_accepts_verified_offline_cache(self):
        self.check_cmake_cache("complete", True)

    def test_libsvm_header_version_is_pinned(self):
        path = ROOT / "third_party" / "libsvm" / "svm.h"
        self.assertTrue(path.is_file(), "missing LIBSVM header")
        self.assertRegex(path.read_text(encoding="utf-8"), r"#define\s+LIBSVM_VERSION\s+337\b")

    def test_relocated_historical_files_preserve_original_bytes(self):
        expected = {
            "examples/regression/train_sequences.tsv": "47450adddd72cf47facc346c583f938dad46e28194b22d12bcfb52f21e818fd4",
            "examples/regression/test_sequences.tsv": "35c46989c7eadb8ecd4c63ae52e085a2dd1f5b7a51adf5c8f8547c41c778db61",
            "legacy/paired_sequence_analyzer.cpp": "a1380f6217d1cdd8fd40f1d69ce7a34f12d1451e6cf0ddfae7a1cad8ebb40318",
        }
        for name, digest in expected.items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(), digest)

    def test_project_version_is_pinned(self):
        content = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertRegex(content, r"project\(StringKernel4LIBSVM\s+VERSION\s+1\.0\.0\b")

    def test_test_framework_archive_is_pinned(self):
        path = ROOT / "cmake" / "Dependencies.cmake"
        self.assertTrue(path.is_file(), "missing test dependency declaration")
        content = path.read_text(encoding="utf-8")
        self.assertIn("https://codeload.github.com/catchorg/Catch2/tar.gz/refs/tags/v3.16.1", content)
        self.assertIn(f"SHA256={CATCH2_SHA256}", content)
        self.assertNotRegex(content, re.compile(r"\bGIT_TAG\s+(main|master)\b"))


if __name__ == "__main__":
    unittest.main()
