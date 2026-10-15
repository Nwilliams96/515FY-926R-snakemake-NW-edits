"""Validate and atomically publish one completed PR2 classifier.

The scientific construction steps remain visible as individual Snakemake
rules. This helper is intentionally limited to the infrastructure operation
that independent Snakemake processes cannot coordinate themselves: locking a
shared destination and replacing it only with a fully validated artifact.
"""

# AI-assisted change — see provenance ledger linked from README.md, commits 8d63a81 and 1a25e17.
# Purpose: prevent concurrent projects from exposing a partial PR2 classifier.
# Human review status: pending.

import fcntl
import os
import shutil
import subprocess
from pathlib import Path


def classifier_is_valid(path, log_handle):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return False
    result = subprocess.run(
        ["qiime", "tools", "validate", str(path)],
        stdout=log_handle,
        stderr=subprocess.STDOUT,
    )
    return result.returncode == 0


def write_marker(marker, classifier):
    marker = Path(marker)
    marker.parent.mkdir(parents=True, exist_ok=True)
    temporary = marker.with_name(marker.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        f"classifier\t{Path(classifier).resolve()}\n"
        f"size_bytes\t{Path(classifier).stat().st_size}\n",
        encoding="utf-8",
    )
    os.replace(temporary, marker)


def publish_classifier(source, destination, lock_path, marker, log_path):
    source = Path(source)
    destination = Path(destination)
    lock_path = Path(lock_path)
    marker = Path(marker)
    log_path = Path(log_path)

    destination.parent.mkdir(parents=True, exist_ok=True)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as log_handle:
        if source.resolve() != destination.resolve() and not classifier_is_valid(
            source, log_handle
        ):
            raise RuntimeError(f"Staged PR2 classifier failed QIIME validation: {source}")

        log_handle.write(f"Waiting for shared PR2 publication lock: {lock_path}\n")
        log_handle.flush()
        with lock_path.open("a+b") as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            log_handle.write("Acquired shared PR2 publication lock\n")
            log_handle.flush()

            if classifier_is_valid(destination, log_handle):
                log_handle.write(f"Reusing validated PR2 classifier: {destination}\n")
                write_marker(marker, destination)
                return

            if source.resolve() == destination.resolve():
                raise RuntimeError(
                    "The configured pre-existing PR2 classifier is invalid: "
                    f"{destination}"
                )

            temporary = destination.with_name(
                destination.name + f".publishing.{os.getpid()}"
            )
            try:
                with source.open("rb") as input_handle, temporary.open("wb") as output:
                    shutil.copyfileobj(input_handle, output)
                    output.flush()
                    os.fsync(output.fileno())
                if not classifier_is_valid(temporary, log_handle):
                    raise RuntimeError(
                        "Copied PR2 classifier failed validation before publication"
                    )
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)

            if not classifier_is_valid(destination, log_handle):
                raise RuntimeError("Published PR2 classifier failed QIIME validation")
            write_marker(marker, destination)
            log_handle.write(f"Published PR2 classifier atomically: {destination}\n")


if "snakemake" in globals():
    publish_classifier(
        snakemake.input.classifier,
        snakemake.params.classifier,
        snakemake.params.lock,
        snakemake.output.marker,
        snakemake.log[0],
    )
