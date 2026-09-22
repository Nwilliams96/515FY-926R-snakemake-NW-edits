"""Create or validate the shared PR2 classifier without cross-run races."""

import fcntl
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path


def checksum(path, algorithm="sha256", chunk_size=1024 * 1024):
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_verified(url, destination, expected_checksum):
    destination = Path(destination)
    with urllib.request.urlopen(url) as response, destination.open("wb") as output:
        shutil.copyfileobj(response, output)
    observed = checksum(destination)
    if observed != expected_checksum:
        raise ValueError(
            "SHA256 checksum mismatch for PR2 reference: "
            f"expected {expected_checksum}, received {observed}"
        )


def reformat_reference(input_path, clean_path, taxonomy_path):
    feature_number = 0
    with Path(input_path).open(encoding="utf-8") as source, Path(clean_path).open(
        "w", encoding="utf-8"
    ) as clean, Path(taxonomy_path).open("w", encoding="utf-8") as taxonomy:
        for line in source:
            if line.startswith(">"):
                feature_number += 1
                feature_id = f"feature_{feature_number}"
                clean.write(f">{feature_id}\n")
                taxonomy.write(f"{feature_id}\t{line[1:].strip()}\n")
            else:
                clean.write(line)
    if feature_number == 0:
        raise ValueError("The downloaded PR2 reference contains no sequences")


def run(command, log_handle):
    log_handle.write("$ " + " ".join(map(str, command)) + "\n")
    log_handle.flush()
    subprocess.run(
        [str(part) for part in command],
        check=True,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
    )


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
    temporary.replace(marker)


def build_classifier(root, params, threads, log_handle):
    source_gz = root / "pr2.fasta.gz"
    source_fasta = root / "pr2.fasta"
    clean_fasta = root / "pr2.clean.fasta"
    taxonomy_tsv = root / "pr2.taxonomy.tsv"
    clean_qza = root / "pr2.clean.qza"
    taxonomy_qza = root / "pr2.taxonomy.qza"
    culled_qza = root / "pr2.clean.culled.qza"
    derep_qza = root / "pr2.clean.culled.derep.qza"
    derep_taxonomy_qza = root / "pr2.taxonomy.derep.qza"
    sliced_qza = root / "pr2.sliced.qza"
    extraction_stats_qza = root / "pr2.read-extraction-stats.qza"
    sliced_derep_qza = root / "pr2.sliced.derep.qza"
    sliced_taxonomy_qza = root / "pr2.sliced.taxonomy.qza"
    classifier_qza = root / "pr2.classifier.qza"

    download_verified(params.source_url, source_gz, params.source_checksum)
    with gzip.open(source_gz, "rb") as compressed, source_fasta.open("wb") as output:
        shutil.copyfileobj(compressed, output)
    reformat_reference(source_fasta, clean_fasta, taxonomy_tsv)

    run(
        [
            "qiime", "tools", "import",
            "--type", "FeatureData[Sequence]",
            "--input-path", clean_fasta,
            "--output-path", clean_qza,
        ],
        log_handle,
    )
    run(
        [
            "qiime", "tools", "import",
            "--type", "FeatureData[Taxonomy]",
            "--input-format", "HeaderlessTSVTaxonomyFormat",
            "--input-path", taxonomy_tsv,
            "--output-path", taxonomy_qza,
        ],
        log_handle,
    )
    run(
        [
            "qiime", "rescript", "cull-seqs",
            "--i-sequences", clean_qza,
            "--o-clean-sequences", culled_qza,
        ],
        log_handle,
    )
    run(
        [
            "qiime", "rescript", "dereplicate",
            "--i-sequences", culled_qza,
            "--i-taxa", taxonomy_qza,
            "--p-mode", "uniq",
            "--p-rank-handles", "disable",
            "--o-dereplicated-sequences", derep_qza,
            "--o-dereplicated-taxa", derep_taxonomy_qza,
        ],
        log_handle,
    )
    run(
        [
            "qiime", "feature-classifier", "extract-reads",
            "--i-sequences", derep_qza,
            "--p-f-primer", params.fwd_primer,
            "--p-r-primer", params.rev_primer,
            "--p-n-jobs", threads,
            "--p-read-orientation", "forward",
            "--o-reads", sliced_qza,
            "--o-read-extraction-stats", extraction_stats_qza,
        ],
        log_handle,
    )
    run(
        [
            "qiime", "rescript", "dereplicate",
            "--i-sequences", sliced_qza,
            "--i-taxa", derep_taxonomy_qza,
            "--p-mode", "uniq",
            "--o-dereplicated-sequences", sliced_derep_qza,
            "--o-dereplicated-taxa", sliced_taxonomy_qza,
        ],
        log_handle,
    )
    run(
        [
            "qiime", "feature-classifier", "fit-classifier-naive-bayes",
            "--i-reference-reads", sliced_derep_qza,
            "--i-reference-taxonomy", sliced_taxonomy_qza,
            "--o-classifier", classifier_qza,
        ],
        log_handle,
    )
    if not classifier_is_valid(classifier_qza, log_handle):
        raise RuntimeError("The newly built PR2 classifier failed QIIME validation")
    return classifier_qza


def ensure_classifier(snakemake_object):
    classifier = Path(str(snakemake_object.params.classifier))
    lock_path = Path(str(snakemake_object.params.lock))
    marker = Path(str(snakemake_object.output.marker))
    log_path = Path(str(snakemake_object.log[0]))
    classifier.parent.mkdir(parents=True, exist_ok=True)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as log_handle, lock_path.open(
        "a+", encoding="utf-8"
    ) as lock_handle:
        log_handle.write(f"Waiting for shared PR2 lock: {lock_path}\n")
        log_handle.flush()
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        log_handle.write("Acquired shared PR2 lock\n")
        log_handle.flush()

        if classifier_is_valid(classifier, log_handle):
            log_handle.write(f"Reusing validated PR2 classifier: {classifier}\n")
            write_marker(marker, classifier)
            return

        if not bool(snakemake_object.params.allow_build):
            raise FileNotFoundError(
                "The configured pre-existing PR2 classifier is missing or invalid: "
                f"{classifier}"
            )

        log_handle.write(
            "No valid PR2 classifier was found; building it in a private staging "
            "directory.\n"
        )
        log_handle.flush()
        with tempfile.TemporaryDirectory(
            prefix=".pr2-build-", dir=classifier.parent
        ) as temporary_directory:
            staged_classifier = build_classifier(
                Path(temporary_directory),
                snakemake_object.params,
                int(snakemake_object.threads),
                log_handle,
            )
            os.replace(staged_classifier, classifier)

        if not classifier_is_valid(classifier, log_handle):
            raise RuntimeError(
                "The atomically published PR2 classifier failed QIIME validation"
            )
        write_marker(marker, classifier)
        log_handle.write(f"Published PR2 classifier atomically: {classifier}\n")


if "snakemake" in globals():
    ensure_classifier(snakemake)
