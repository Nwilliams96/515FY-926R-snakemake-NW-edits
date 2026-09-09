"""Download a reference artifact atomically and verify its checksum."""

import hashlib
import shutil
import urllib.request
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


def download_verified(url, destination, expected_checksum, algorithm="md5"):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.is_file() and checksum(destination, algorithm) == expected_checksum:
        return

    partial = destination.with_name(destination.name + ".part")
    partial.unlink(missing_ok=True)
    try:
        with urllib.request.urlopen(url) as response, open(partial, "wb") as output:
            shutil.copyfileobj(response, output)
        observed_checksum = checksum(partial, algorithm)
        if observed_checksum != expected_checksum:
            raise ValueError(
                f"{algorithm.upper()} checksum mismatch for {url}: "
                f"expected {expected_checksum}, received {observed_checksum}"
            )
        partial.replace(destination)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise


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
