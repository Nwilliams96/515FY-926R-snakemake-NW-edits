"""Resolve classifier choices and paths; scientific steps stay in Snakemake.

"""

import hashlib
import json
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path
from urllib.parse import urlparse


SILVA_REGION_NAME = "SILVA_144_SSURef_NR99_uniform_classifier_V4V5-515f-926r.qza"
SILVA_FULL_NAME = (
    "SILVA_144_SSURef_NR99_full-length_prokaryotes_"
    "rank-propagation_uniq_genus_uniform_classifier.qza"
)
SILVA_BASE = "https://www.arb-silva.de/fileadmin/silva_databases/release_144/QIIME2/2026.7/SSU/"
SILVA_PRESETS = {
    SILVA_BASE + "V4V5-515f-926r/uniform/" + SILVA_REGION_NAME:
        (SILVA_REGION_NAME, "f7757b01eb82e0ac78bf06427e410095"),
    SILVA_BASE + "full-length/uniform/" + SILVA_FULL_NAME:
        (SILVA_FULL_NAME, "e16b9bc78dafa2f344563e570bbf25ba"),
}


def resolve_qiime_environment(specification, source_path):
    """Resolve YAML or an installed name without creating or changing an env."""
    specification = str(specification).strip()
    if specification.endswith((".yaml", ".yml")):
        try:
            content = Path(source_path(specification)).read_bytes()
        except OSError as error:
            raise ValueError(f"Cannot read QIIME 2 environment YAML: {specification}") from error
        # Preserve existing YAML-based classifier cache names exactly.
        return dict(kind="yaml", requested=specification, specification=specification,
                    identity=hashlib.sha256(content).hexdigest()[:12])

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", specification):
        raise ValueError("qiime2version must be a .yaml/.yml file or an installed Conda environment name")
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if not conda:
        raise ValueError("Conda is unavailable; provide a QIIME 2 YAML file or make Conda available to resolve the environment name")

    def conda_json(*arguments):
        try:
            result = subprocess.run([conda, *arguments], check=True, capture_output=True,
                                    text=True, timeout=30)
            return json.loads(result.stdout)
        except (OSError, subprocess.SubprocessError, ValueError) as error:
            raise ValueError(f"Cannot inspect Conda environment {specification!r}; check Conda availability and permissions") from error

    info = conda_json("info", "--json")
    if not isinstance(info, dict):
        raise ValueError("Conda returned an invalid environment inventory")
    prefixes = info.get("envs", [])
    if not isinstance(prefixes, list):
        raise ValueError("Conda returned an invalid environment list")
    candidates = {str(Path(prefix).resolve()) for prefix in prefixes
                  if isinstance(prefix, str) and Path(prefix).name == specification}
    if specification == "base" and info.get("root_prefix"):
        candidates = {str(Path(info["root_prefix"]).resolve())}
    if len(candidates) != 1:
        reason = "was not found" if not candidates else "is ambiguous across Conda prefixes"
        raise ValueError(f"QIIME 2 Conda environment {specification!r} {reason}; use a YAML definition instead")
    prefix = candidates.pop()
    if not Path(prefix, "conda-meta").is_dir():
        raise ValueError(f"QIIME 2 environment is not an installed Conda prefix: {prefix}")
    packages = conda_json("list", "--prefix", prefix, "--json")
    if not isinstance(packages, list) or not packages or any(
        not isinstance(package, dict) or not package.get("name") or not package.get("version")
        for package in packages
    ):
        raise ValueError(f"Cannot fingerprint installed packages in QIIME 2 environment {specification!r}")
    # Normalize ordering and omit machine-specific paths from the fingerprint.
    # Include Conda and pip inventory entries; a name alone is not a recipe.
    inventory = sorted([
        {key: package.get(key, "") for key in ("name", "version", "build_string", "channel", "platform")}
        for package in packages
    ], key=lambda package: json.dumps(package, sort_keys=True))
    identity = "named-" + hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest()[:12]
    versions = {package["name"]: package["version"] for package in packages
                if package["name"] in {"qiime2", "scikit-learn", "q2-dada2", "q2-feature-classifier", "rescript"}}
    # Activate the exact prefix inspected above, avoiding later name ambiguity.
    return dict(kind="named", requested=specification, specification=prefix,
                identity=identity, software_versions=versions)


def reusable_classifier(path, expected_md5=None):
    """Inspect type/ZIP integrity at DAG time; QIIME validates in a rule."""
    try:
        with zipfile.ZipFile(path) as archive:
            # Provenance also contains metadata.yaml files. Only inspect the
            # archive's root metadata, not metadata inside its provenance DAG.
            metadata = [n for n in archive.namelist()
                        if len(n.split("/")) == 2 and n.endswith("/metadata.yaml")]
            if len(metadata) != 1 or archive.testzip() is not None:
                return False
            text = archive.read(metadata[0]).decode("utf-8")
            if not re.search(r"^type:\s*TaxonomicClassifier\s*$", text, re.M):
                return False
        if expected_md5:
            digest = hashlib.md5()
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            if digest.hexdigest().lower() != expected_md5.lower():
                return False
        return True
    except (OSError, ValueError, zipfile.BadZipFile, UnicodeError):
        return False


def resolve_database(config, family):
    """Availability wins; otherwise a nonblank URL selects download."""
    lower = family.lower()
    version_key = "SILVAversion" if family == "SILVA" else "PR2version"
    version = str(config.get(version_key, "144" if family == "SILVA" else "5.1.1"))
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)*", version):
        raise ValueError(f"{version_key} must be a release number such as 144 or 5.1.1")
    root = os.path.normpath(os.path.expanduser(config["database_dir"]))
    directory = os.path.normpath(os.path.expanduser(
        config.get(lower + "_database_dir") or os.path.join(root, "classification", family)
    ))
    existing = bool(config.get("use_preexisting_" + lower + "_database", False))
    url = str(config.get(family + "_download_link") or "").strip()
    mode = "existing" if existing else ("download" if url else "train")
    region = config.get(lower + "_reference_mode", "amplicon")
    if region not in ("amplicon", "full_length"):
        raise ValueError(f"{lower}_reference_mode must be amplicon or full_length")
    primers = config["fwdPrimer"] + "_" + config["revPrimer"]
    default_primers = primers == "GTGYCAGCMGCCGCGGTAA_CCGYCAATTYMTTTRAGTTT"
    if family == "SILVA":
        legacy = []
        if version == "144":
            legacy.append(SILVA_REGION_NAME if default_primers else SILVA_FULL_NAME)
            if not default_primers:
                # Retain read-only reuse of the old local spelling. Downloads
                # use the verified official filename with an underscore.
                legacy.append(SILVA_FULL_NAME.replace("prokaryotes_", "prokaryotes-"))
        legacy.append("silva-ssu-nr99-tax-dereplicated-sliced_" + primers +
                      "_dereplicated_final_classifier_USE_ME.qza")
    else:
        stem = "pr2_version_" + version + "_SSU_dada2.clean.culled.derep-sliced_" + primers
        legacy = [stem + "_dereplicated_final_classifier_qiime2-2026.7.qza",
                  stem + "_dereplicated_final_classifier_USE_ME.qza"]
    checksum = None
    if mode == "existing":
        supplied = str(config.get(lower + "_classifier") or "").strip()
        candidates = [os.path.join(directory, name) for name in legacy]
        classifier = os.path.expanduser(supplied) if supplied else next(
            (p for p in candidates if os.path.isfile(p)), candidates[0]
        )
        if not os.path.isfile(classifier):
            raise ValueError(f"Pre-existing {family} classifier was not found: {classifier}")
        if not reusable_classifier(classifier):
            raise ValueError(f"Pre-existing {family} file is not an intact TaxonomicClassifier: {classifier}")
        directory = os.path.dirname(os.path.abspath(classifier))
    elif mode == "download":
        if urlparse(url).scheme not in ("http", "https") or not urlparse(url).netloc:
            raise ValueError(f"{family}_download_link must be an HTTP(S) URL")
        if family == "SILVA" and url in SILVA_PRESETS and version != "144":
            raise ValueError("The official SILVA 144 preset requires SILVAversion: '144'")
        checksum = str(config.get(family + "_download_md5") or "").strip() or None
        if checksum and not re.fullmatch(r"[a-fA-F0-9]{32}", checksum):
            raise ValueError(f"{family}_download_md5 must contain 32 hexadecimal characters")
        if checksum:
            checksum = checksum.lower()
        # Keep the current official filenames, but do not infer a checksum
        # for arbitrary URLs or infer the release of a custom classifier.
        if family == "SILVA" and version == "144" and url in SILVA_PRESETS:
            name, published = SILVA_PRESETS[url]
            if checksum and checksum.lower() != published:
                raise ValueError("The supplied checksum disagrees with the official SILVA preset")
            checksum = published
        else:
            identity = hashlib.sha256((url + "\n" + (checksum or "")).encode()).hexdigest()[:16]
            name = f"{lower}-{version}-download-{identity}.qza"
        classifier = os.path.join(directory, name)
    else:
        # Include recipe and environment definition identity to separate
        # training from URL downloads and from earlier PR2 preparation recipes.
        environment = config.get("_database_environment_id", "qiime2-2026.7")
        scope = region + ("-" + primers if region == "amplicon" else "")
        recipe = "rescript-v1-all-domains-min1200-900-1400" if family == "SILVA" else "rescript-v1-nine-ranks"
        name = f"{lower}-{version}-{recipe}-{scope}-{environment}.qza"
        classifier = os.path.join(directory, name)
    identity = hashlib.sha256(os.path.abspath(classifier).encode()).hexdigest()[:16]
    private = os.path.join(".snakemake", "database-build", family, identity)
    ready = os.path.join(".snakemake", "database-ready", lower + "-" + identity + ".ready")
    classifier_region = "as supplied; verify compatibility with your primers"
    if mode == "download" and family == "SILVA" and version == "144" and url in SILVA_PRESETS:
        classifier_region = "amplicon (515Y/926R)" if SILVA_PRESETS[url][0] == SILVA_REGION_NAME else "full_length"
    # Record the effective, visible Snakemake recipe rather than making the
    # report infer scientific settings from a long cache filename. These
    # settings describe local training only, never an externally supplied model.
    training_settings = {}
    if mode == "train":
        training_settings = {
            "Acquisition": "RESCRIPt get-silva-data, SSURef_NR99" if family == "SILVA" else "RESCRIPt get-pr2-data",
            "Reference scope": "Bacteria, Archaea, Eukaryota" if family == "SILVA" else "All nine PR2 ranks; :plas labels retained",
            "Culling": "Remove sequences with at least 5 degenerate bases or homopolymers of at least 8 bases",
            "Dereplication": "uniq, before training and again after primer extraction when enabled",
            "Primer extraction": "Configured forward/reverse primers; forward orientation" if region == "amplicon" else "Skipped (full-length training)",
            "Training": "QIIME feature-classifier fit-classifier-naive-bayes",
            "Environment definition ID": config.get("_database_environment_id", "qiime2-2026.7"),
        }
        if family == "SILVA":
            training_settings["Minimum lengths (nt)"] = "Bacteria 1200; Archaea 900; Eukaryota 1400"
        if region == "amplicon":
            training_settings["Forward primer"] = config["fwdPrimer"]
            training_settings["Reverse primer"] = config["revPrimer"]
    return dict(mode=mode, version=version, classifier=classifier, directory=directory,
                url=url if mode == "download" else "", checksum=checksum,
                reference_mode=region if mode == "train" else classifier_region,
                training_settings=training_settings,
                build_root=private, staged=os.path.join(private, "classifier.qza"),
                ready=ready, lock=(ready + ".lock" if mode == "existing" else classifier + ".lock"))
