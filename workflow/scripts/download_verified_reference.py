"""Download a shared reference safely and verify its checksum."""

import fcntl
import hashlib
import os
import shutil
import tempfile
import urllib.request
from contextlib import contextmanager
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


def download_verified(url, destination, expected_checksum, algorithm="md5"):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    with destination_lock(destination):
        # A different pipeline may have completed the download while this run
        # waited for the lock. Revalidate only after acquiring it.
        if (
            destination.is_file()
            and checksum(destination, algorithm) == expected_checksum
        ):
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
                with urllib.request.urlopen(url) as response:
                    shutil.copyfileobj(response, output)
                output.flush()
                os.fsync(output.fileno())

            observed_checksum = checksum(partial, algorithm)
            if observed_checksum != expected_checksum:
                raise ValueError(
                    f"{algorithm.upper()} checksum mismatch for {url}: "
                    f"expected {expected_checksum}, received {observed_checksum}"
                )

            # The old valid reference remains readable until the fully
            # verified replacement is atomically published.
            os.replace(partial, destination)
            partial = None
        finally:
            if partial is not None:
                partial.unlink(missing_ok=True)


if "snakemake" in globals():
    checksum_algorithm = getattr(snakemake.params, "checksum_algorithm", "md5")
    expected_checksum = getattr(snakemake.params, "checksum", None)
    if expected_checksum is None:
        expected_checksum = snakemake.params["md5"]
    download_verified(
        snakemake.params["url"],
        snakemake.output[0],
        expected_checksum,
        checksum_algorithm,
    )
