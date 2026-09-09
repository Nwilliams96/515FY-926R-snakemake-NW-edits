"""Download a reference artifact atomically and verify its published MD5."""

import hashlib
import shutil
import urllib.request
from pathlib import Path


def md5sum(path, chunk_size=1024 * 1024):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_verified(url, destination, expected_md5):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.is_file() and md5sum(destination) == expected_md5:
        return

    partial = destination.with_name(destination.name + ".part")
    partial.unlink(missing_ok=True)
    try:
        with urllib.request.urlopen(url) as response, open(partial, "wb") as output:
            shutil.copyfileobj(response, output)
        observed_md5 = md5sum(partial)
        if observed_md5 != expected_md5:
            raise ValueError(
                f"Checksum mismatch for {url}: expected {expected_md5}, "
                f"received {observed_md5}"
            )
        partial.replace(destination)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise


if "snakemake" in globals():
    download_verified(
        snakemake.params["url"],
        snakemake.output[0],
        snakemake.params["md5"],
    )
