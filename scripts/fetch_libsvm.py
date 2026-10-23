#!/usr/bin/env python3
"""Prepare the exact official LIBSVM sources without committing a solver copy."""

from __future__ import annotations

import argparse
import hashlib
from http.client import HTTPException
import json
from pathlib import Path
import re
import sys
import tempfile
from urllib import request


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = ROOT / "cmake" / "libsvm.json"
DEFAULT_SOURCE_DIR = ROOT / "third_party" / "libsvm"
REQUIRED_FILES = {"svm.cpp", "svm.h", "svm-train.c", "svm-predict.c", "README", "COPYRIGHT"}


class BootstrapError(RuntimeError):
    """The dependency could not be prepared and must not be used."""


def load_manifest(path: Path) -> dict:
    """Reject incomplete pins, unsafe names, and floating or unofficial URLs."""
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(manifest, dict):
            raise ValueError("expected an object")
        for key in ("version", "tag"):
            if not isinstance(manifest.get(key), str) or not manifest[key].strip():
                raise ValueError(f"missing {key}")
        commit = manifest.get("commit")
        if not isinstance(commit, str) or not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise ValueError("commit must be a complete lowercase Git commit hash")
        files = manifest.get("files")
        if not isinstance(files, dict) or set(files) != REQUIRED_FILES:
            raise ValueError("expected the complete six-file LIBSVM source set")
        for name, entry in files.items():
            if not isinstance(entry, dict):
                raise ValueError(f"invalid entry for {name}")
            digest = entry.get("sha256")
            if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
                raise ValueError(f"invalid SHA-256 for {name}")
            expected_url = f"https://raw.githubusercontent.com/cjlin1/libsvm/{commit}/{name}"
            if entry.get("url") != expected_url:
                raise ValueError(f"expected the pinned official source URL for {name}")
        return manifest
    except (OSError, UnicodeError, ValueError) as exc:
        raise BootstrapError(f"Invalid LIBSVM manifest {path}: {exc}") from exc


def verify_bytes(content: bytes, digest: str, name: str) -> None:
    if hashlib.sha256(content).hexdigest() != digest:
        raise BootstrapError(f"SHA-256 mismatch for {name}; remove the modified file and "
                             "rerun python3 scripts/fetch_libsvm.py.")


def ensure_libsvm(source_dir: Path, manifest_path: Path = DEFAULT_MANIFEST) -> None:
    """Verify a cache offline, or stage and verify all missing files before use."""
    manifest = load_manifest(Path(manifest_path))
    source_dir = Path(source_dir)
    if source_dir.exists() and not source_dir.is_dir():
        raise BootstrapError(f"Invalid LIBSVM source directory {source_dir}; remove the file "
                             "or choose a directory with --source-dir.")
    missing = []
    for name, entry in manifest["files"].items():
        path = source_dir / name
        if path.exists():
            if not path.is_file():
                raise BootstrapError(f"Invalid LIBSVM file {path}; remove it and rerun the bootstrap.")
            verify_bytes(path.read_bytes(), entry["sha256"], name)
        else:
            missing.append(name)
    if not missing:
        return

    source_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".libsvm-fetch-", dir=source_dir.parent) as directory:
        staging = Path(directory)
        for name in missing:
            entry = manifest["files"][name]
            try:
                with request.urlopen(entry["url"], timeout=30) as response:
                    content = response.read()
            except (OSError, HTTPException) as exc:
                raise BootstrapError(f"Could not download {name}: {exc}. Check network access "
                                     "and rerun python3 scripts/fetch_libsvm.py.") from exc
            verify_bytes(content, entry["sha256"], name)
            (staging / name).write_bytes(content)
        source_dir.mkdir(parents=True, exist_ok=True)
        for name in missing:
            (staging / name).replace(source_dir / name)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR,
                        help="cache destination (default: third_party/libsvm)")
    options = parser.parse_args(argv)
    try:
        ensure_libsvm(options.source_dir)
    except (BootstrapError, OSError) as exc:
        print(f"LIBSVM bootstrap failed: {exc}", file=sys.stderr)
        return 1
    print(f"Verified LIBSVM sources: {options.source_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
