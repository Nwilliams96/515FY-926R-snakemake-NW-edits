"""Download a shared reference safely and verify its checksum."""


import fcntl
import hashlib
import os
import shutil
import re
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from contextlib import contextmanager, nullcontext
from pathlib import Path


def checksum(path, algorithm="md5", chunk_size=1024 * 1024):
    digest = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5sum(path, chunk_size=1024 * 1024):
    """Retain the original helper for callers using published MD5 values."""
    return checksum(path, "md5", chunk_size)


@contextmanager
def destination_lock(destination):
    """Serialize independent pipeline runs publishing the same reference."""
    lock_path = destination.with_name(destination.name + ".lock")
    with open(lock_path, "a+b") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_handle, fcntl.LOCK_UN)


def validate_classifier_file(path, log_handle):
    """Reject HTML, wrong artifact types, and corrupt QIIME archives."""
    with zipfile.ZipFile(path) as archive:
        roots = [name for name in archive.namelist()
                 if len(name.split("/")) == 2 and name.endswith("/metadata.yaml")]
        if len(roots) != 1 or not re.search(
            r"^type:\s*TaxonomicClassifier\s*$",
            archive.read(roots[0]).decode("utf-8"), re.M
        ):
            raise ValueError("The download is not a QIIME TaxonomicClassifier artifact")
    subprocess.run(["qiime", "tools", "validate", str(path)], check=True,
                   stdout=log_handle, stderr=subprocess.STDOUT)


def download_verified(url, destination, expected_checksum=None, algorithm="md5",
                      validate_classifier=False, log_path=None):
    # Custom pre-trained classifier URLs may lack a published checksum, but
    # must pass QIIME type and integrity validation. Raw downloads retain the
    # original mandatory-checksum contract (including the BBsplit references).
    if not expected_checksum and not validate_classifier:
        raise ValueError("Raw reference downloads require a published checksum")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    log_context = open(log_path, "a", encoding="utf-8") if log_path else nullcontext(sys.stderr)

    with log_context as log_handle, destination_lock(destination):
        # A different pipeline may have completed the download while this run
        # waited for the lock. Revalidate only after acquiring it.
        if (
            destination.is_file()
            and (not expected_checksum or checksum(destination, algorithm).lower() == expected_checksum.lower())
        ):
            try:
                if validate_classifier:
                    validate_classifier_file(destination, log_handle)
            except (ValueError, OSError, UnicodeError, zipfile.BadZipFile, subprocess.CalledProcessError):
                # An invalid cached file is never mistaken for completion.
                pass
            else:
                log_handle.write(f"Reusing validated reference: {destination}\n")
                return

        partial = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=destination.parent,
                prefix=f".{destination.name}.",
                suffix=".part",
                delete=False,
            ) as output:
                partial = Path(output.name)
                with urllib.request.urlopen(url, timeout=120) as response:
                    shutil.copyfileobj(response, output)
                output.flush()
                os.fsync(output.fileno())

            if expected_checksum:
                observed_checksum = checksum(partial, algorithm)
                if observed_checksum.lower() != expected_checksum.lower():
                    raise ValueError(
                        f"{algorithm.upper()} checksum mismatch for {url}: "
                        f"expected {expected_checksum}, received {observed_checksum}"
                    )
            else:
                log_handle.write("No publisher checksum supplied; validating classifier type and archive.\n")
            if validate_classifier:
                validate_classifier_file(partial, log_handle)

            # The old valid reference remains readable until the fully
            # verified replacement is atomically published.
            os.replace(partial, destination)
            log_handle.write(f"Published validated reference: {destination}\n")
            partial = None
        finally:
            if partial is not None:
                partial.unlink(missing_ok=True)


if "snakemake" in globals():
    checksum_algorithm = getattr(snakemake.params, "checksum_algorithm", "md5")
    expected_checksum = getattr(snakemake.params, "checksum", None)
    if expected_checksum is None:
        expected_checksum = getattr(snakemake.params, "md5", None)
    download_verified(
        snakemake.params["url"],
        snakemake.output[0],
        expected_checksum,
        checksum_algorithm,
        validate_classifier=getattr(snakemake.params, "validate_classifier", False),
        log_path=str(snakemake.log[0]) if snakemake.log else None,
    )
