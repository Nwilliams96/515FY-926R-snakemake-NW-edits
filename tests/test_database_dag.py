"""Isolated dry-runs of the real workflow; no reference downloads or QIIME runs.

The workflow copy receives synthetic classifiers and a minimal generated-config
fixture. Tests are skipped when Snakemake is not installed in this interpreter.
"""

import importlib.util
import itertools
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is not installed")
class DatabaseDagTests(unittest.TestCase):
    def setup_project(self, root):
        root = Path(root)
        shutil.copytree(ROOT / "workflow", root / "workflow")
        config_dir = root / "config"
        (config_dir / "schemas").mkdir(parents=True)
        # Configuration-builder tests exercise its full public schema. Here
        # only valid sample identifiers are needed to parse the actual DAG.
        schema_version = '$schema: "https://json-schema.org/draft/2020-12/schema"\n'
        (config_dir / "schemas/config.schema.yml").write_text(schema_version + "type: object\n")
        (config_dir / "schemas/samples.schema.yml").write_text(
            schema_version + "type: object\nproperties:\n  sample:\n    type: string\nrequired: [sample]\n"
        )
        (config_dir / "samples.tsv").write_text("sample\n001\n")
        (config_dir / "prok_and_euk_SSU_amplicon_molarities.tsv").write_text(
            "sample_type\tamount_pM\n16S\t1000\n18S\t1000\n"
        )
        shared = root / "shared databases"
        (shared / "bbsplit-db/EUK-PROK-bbsplit-db").mkdir(parents=True)
        return dict(
            samplesheet="config/samples.tsv", studyName="DAG-fixture",
            projectName="DAG-fixture", database_dir=str(shared),
            rawdatadir=str(root / "raw"), R1file_ending=".R1.fastq.gz",
            R2file_ending=".R2.fastq.gz", use_internal_standards=False,
            use_preexisting_bbsplit_database=True,
            trunclens=dict(truncR1=220, truncR2=180),
            fwdPrimer="GTGYCAGCMGCCGCGGTAA", revPrimer="CCGYCAATTYMTTTRAGTTT",
            qiime2version="../envs/rachis-qiime2-linux-64-2026.7.yml",
            SILVAversion="138.2", PR2version="5.1.1",
        )

    def set_mode(self, root, config, family, mode):
        lower = family.lower()
        config["use_preexisting_" + lower + "_database"] = mode == "existing"
        if mode == "existing":
            artifact = Path(root) / "provided" / (family + ".qza")
            artifact.parent.mkdir(exist_ok=True)
            with zipfile.ZipFile(artifact, "w") as archive:
                archive.writestr("fixture/metadata.yaml", "type: TaxonomicClassifier\n")
                archive.writestr("fixture/data/classifier.pkl", b"fixture")
                archive.writestr(
                    "fixture/provenance/artifacts/ancestor/metadata.yaml",
                    "type: FeatureData[Sequence]\n",
                )
            config[lower + "_classifier"] = str(artifact)
        elif mode == "download":
            config[family + "_download_link"] = f"https://references.example/{family}.qza"

    def dry_run(self, root, config, targets=("ensure_silva_classifier", "ensure_pr2_classifier")):
        (Path(root) / "config/config.yml").write_text(json.dumps(config))
        process = subprocess.run(
            [sys.executable, "-m", "snakemake", "--snakefile", "workflow/Snakefile",
             "--dry-run", "--forceall", "--nolock", "--shared-fs-usage", "input-output",
             "persistence", "--cores", "8", "--printshellcmds", *targets],
            cwd=root, text=True, capture_output=True, timeout=60,
        )
        output = process.stdout + process.stderr
        self.assertEqual(process.returncode, 0, output)
        return output

    def test_all_nine_choices_schedule_only_the_selected_preparation_route(self):
        with tempfile.TemporaryDirectory() as tmp:
            for silva_mode, pr2_mode in itertools.product(("existing", "download", "train"), repeat=2):
                with self.subTest(SILVA=silva_mode, PR2=pr2_mode):
                    root = Path(tmp) / (silva_mode + "-" + pr2_mode)
                    config = self.setup_project(root)
                    self.set_mode(root, config, "SILVA", silva_mode)
                    self.set_mode(root, config, "PR2", pr2_mode)
                    output = self.dry_run(root, config)
                    for family, mode in (("SILVA", silva_mode), ("PR2", pr2_mode)):
                        download = "rule download_" + family + "_classifier:"
                        train = "rule train_classifier:" if family == "SILVA" else "rule train_classifier_pr2:"
                        extract = "rule extract_primers:" if family == "SILVA" else "rule extract_primers_pr2:"
                        self.assertEqual(download in output, mode == "download", output)
                        self.assertEqual(train in output, mode == "train", output)
                        self.assertEqual(extract in output, mode == "train", output)

    def test_full_length_training_bypasses_primer_extraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.setup_project(tmp)
            config.update(silva_reference_mode="full_length", pr2_reference_mode="full_length")
            output = self.dry_run(tmp, config)
            self.assertIn("rule train_classifier:", output)
            self.assertIn("rule train_classifier_pr2:", output)
            self.assertNotIn("rule extract_primers:", output)
            self.assertNotIn("rule extract_primers_pr2:", output)
            self.assertNotIn("rule dereplicated_sliced_data:", output)
            self.assertNotIn("rule dereplicate_extracted_pr2_reads:", output)

    def test_native_version_capability_check_precedes_raw_read_processing(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.setup_project(tmp)
            config["SILVAversion"] = "144"
            self.set_mode(tmp, config, "PR2", "existing")
            raw = Path(config["rawdatadir"])
            raw.mkdir()
            for direction in ("R1", "R2"):
                (raw / ("001." + direction + ".fastq.gz")).touch()
            output = self.dry_run(tmp, config, ("results/00-trimmed/001.qc.txt",))
            self.assertIn("rule check_silva_reference_capability:", output)
            self.assertIn("rule cutadapt:", output)
            cutadapt = output.split("rule cutadapt:", 1)[1]
            self.assertRegex(cutadapt.split("output:", 1)[0],
                             r"database-build/SILVA/[a-f0-9]+/capability.json")

    def test_completed_compatible_downloads_and_training_are_reused(self):
        spec = importlib.util.spec_from_file_location(
            "database_config", ROOT / "workflow/scripts/database_config.py"
        )
        resolver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(resolver)
        for mode in ("download", "train"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                config = self.setup_project(tmp)
                for family in ("SILVA", "PR2"):
                    self.set_mode(tmp, config, family, mode)
                env = Path(tmp) / "workflow/envs/rachis-qiime2-linux-64-2026.7.yml"
                settings = dict(config, _database_environment_id=hashlib.sha256(env.read_bytes()).hexdigest()[:12])
                for family in ("SILVA", "PR2"):
                    path = Path(resolver.resolve_database(settings, family)["classifier"])
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with zipfile.ZipFile(path, "w") as archive:
                        archive.writestr("fixture/metadata.yaml", "type: TaxonomicClassifier\n")
                        archive.writestr("fixture/data/classifier.pkl", b"completed fixture")
                output = self.dry_run(tmp, config)
                self.assertIn("rule ensure_silva_classifier:", output)
                self.assertIn("rule ensure_pr2_classifier:", output)
                for forbidden in ("download_SILVA_classifier", "download_PR2_classifier",
                                  "train_classifier", "train_classifier_pr2",
                                  "check_silva_reference_capability", "check_pr2_reference_capability"):
                    self.assertNotIn("rule " + forbidden + ":", output)

    def test_full_results_dag_includes_normalized_pr2_inputs_and_matching_species_release(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.setup_project(tmp)
            for family in ("SILVA", "PR2"):
                self.set_mode(tmp, config, family, "existing")
            raw = Path(config["rawdatadir"])
            raw.mkdir()
            for direction in ("R1", "R2"):
                (raw / ("001." + direction + ".fastq.gz")).touch()
            output = self.dry_run(tmp, config, ("all",))
            for rule in ("normalize_canonical_16S_PR2_rank_prefixes",
                         "normalize_unmerged_16S_PR2_rank_prefixes",
                         "normalize_18S_PR2_rank_prefixes",
                         "import_normalized_18S_PR2_taxonomy",
                         "resolve_PR2_plastids", "split_resolved_prok_categories",
                         "generate_pipeline_report", "export_results_for_download"):
                self.assertIn("rule " + rule + ":", output)
            self.assertNotIn("rule download_SILVA_144_species_reference:", output)
            self.assertNotIn("rule assign_SILVA_144_species:", output)
            resolver = output.split("rule resolve_PR2_plastids:", 1)[1].split("output:", 1)[0]
            self.assertIn(".snakemake/classification-exports/canonical-16S/PR2/taxonomy.tsv", resolver)
            supplement = output.split("rule export_unmerged_16S_supplement:", 1)[1].split("output:", 1)[0]
            self.assertIn(".snakemake/classification-exports/unmerged-16S/PR2/taxonomy.tsv", supplement)
            self.assertIn(
                ".snakemake/classification-exports/unmerged-16S/sequences/linked-dna-sequences.fasta",
                supplement,
            )
            self.assertNotIn("/sequences/dna-sequences.fasta", supplement)
            euk_export = output.split("rule export_tax_convert_biom_euk:", 1)[1].split("output:", 1)[0]
            self.assertIn("DAG-fixture_PR2.normalized.classified.qza", euk_export)
            self.assertIn("rule assign_SILVA_species:", output)


if __name__ == "__main__":
    unittest.main()
