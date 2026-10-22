import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("database_config", ROOT / "workflow/scripts/database_config.py")
DATABASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DATABASE)


class QiimeEnvironmentTests(unittest.TestCase):
    def packages(self, version="1.7.1"):
        return [dict(name="qiime2", version="2026.7", build_string="py_0", channel="qiime2"),
                dict(name="scikit-learn", version=version, build_string="py_0", channel="conda-forge")]

    def resolve_named(self, prefix, packages=None, name=None, info=None):
        results = [json.dumps(info or {"envs": [str(prefix)]}),
                   json.dumps(self.packages() if packages is None else packages)]
        with mock.patch.dict(DATABASE.os.environ, {"CONDA_EXE": "/fixture/conda"}), \
             mock.patch.object(DATABASE.subprocess, "run", side_effect=[
                 mock.Mock(stdout=value) for value in results
             ]) as run:
            record = DATABASE.resolve_qiime_environment(name or prefix.name, lambda value: value)
        return record, run

    def test_yaml_preserves_previous_cache_identity_without_running_conda(self):
        with tempfile.TemporaryDirectory() as tmp:
            for ending in (".yml", ".yaml"):
                path = Path(tmp) / ("custom environment" + ending)
                path.write_bytes(b"dependencies:\n- qiime2=2026.7\n")
                with mock.patch.object(DATABASE.subprocess, "run") as run:
                    result = DATABASE.resolve_qiime_environment(path.name, lambda value: Path(tmp) / value)
                self.assertEqual(result["specification"], path.name)
                self.assertEqual(result["kind"], "yaml")
                self.assertEqual(result["identity"], hashlib.sha256(path.read_bytes()).hexdigest()[:12])
                run.assert_not_called()

    def test_missing_yaml_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "Cannot read QIIME 2 environment YAML"):
            DATABASE.resolve_qiime_environment("absent.yml", lambda value: "/missing/" + value)

    def test_named_environment_selects_the_exact_inspected_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "qiime-2026.7"
            (prefix / "conda-meta").mkdir(parents=True)
            record, run = self.resolve_named(prefix)
            self.assertEqual(record["kind"], "named")
            self.assertEqual(record["requested"], prefix.name)
            self.assertEqual(record["specification"], str(prefix.resolve()))
            self.assertEqual(record["software_versions"]["scikit-learn"], "1.7.1")
            self.assertRegex(record["identity"], r"^named-[a-f0-9]{12}$")
            self.assertEqual(run.call_args_list[0].args[0], ["/fixture/conda", "info", "--json"])
            self.assertEqual(run.call_args_list[1].args[0], ["/fixture/conda", "list", "--prefix", str(prefix.resolve()), "--json"])
            self.assertTrue(all(call.kwargs["timeout"] == 30 for call in run.call_args_list))

    def test_package_order_is_irrelevant_but_installed_versions_change_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "qiime"
            (prefix / "conda-meta").mkdir(parents=True)
            original, _ = self.resolve_named(prefix)
            reordered, _ = self.resolve_named(prefix, list(reversed(self.packages())))
            changed, _ = self.resolve_named(prefix, self.packages("1.7.2"))
            self.assertEqual(original["identity"], reordered["identity"])
            self.assertNotEqual(original["identity"], changed["identity"])
            config = dict(database_dir=tmp, fwdPrimer="ACGT", revPrimer="ACGT")
            old = DATABASE.resolve_database(dict(config, _database_environment_id=original["identity"]), "PR2")
            new = DATABASE.resolve_database(dict(config, _database_environment_id=changed["identity"]), "PR2")
            self.assertNotEqual(old["classifier"], new["classifier"])

    def test_base_environment_uses_root_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "miniforge"
            (prefix / "conda-meta").mkdir(parents=True)
            result, _ = self.resolve_named(prefix, name="base", info={"envs": [str(prefix)], "root_prefix": str(prefix)})
            self.assertEqual(result["specification"], str(prefix.resolve()))

    def test_missing_ambiguous_and_invalid_prefixes_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "qiime"
            for info, message in (
                ({"envs": []}, "was not found"),
                ({"envs": [str(prefix), str(Path(tmp) / "other/qiime")]}, "is ambiguous"),
                ({"envs": [str(prefix)]}, "not an installed Conda prefix"),
            ):
                with self.subTest(info=info), self.assertRaisesRegex(ValueError, message):
                    self.resolve_named(prefix, info=info)

    def test_empty_package_inventory_fails_instead_of_using_a_name_only_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "qiime"
            (prefix / "conda-meta").mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "Cannot fingerprint installed packages"):
                self.resolve_named(prefix, [])

    def test_missing_conda_and_failed_inspection_are_clear(self):
        with mock.patch.dict(DATABASE.os.environ, {"CONDA_EXE": ""}), \
             mock.patch.object(DATABASE.shutil, "which", return_value=None), \
             self.assertRaisesRegex(ValueError, "Conda is unavailable"):
            DATABASE.resolve_qiime_environment("qiime", str)
        for error in (subprocess.TimeoutExpired("conda", 30), OSError("unavailable")):
            with mock.patch.dict(DATABASE.os.environ, {"CONDA_EXE": "/fixture/conda"}), \
                 mock.patch.object(DATABASE.subprocess, "run", side_effect=error), \
                 self.assertRaisesRegex(ValueError, "Cannot inspect Conda environment"):
                DATABASE.resolve_qiime_environment("qiime", str)

    def test_unsafe_or_empty_names_are_rejected_before_any_subprocess(self):
        for value in ("", "../qiime", "qiime;echo bad", "qiime\nother", "--name=bad"):
            with self.subTest(value=value), mock.patch.object(DATABASE.subprocess, "run") as run:
                with self.assertRaises(ValueError):
                    DATABASE.resolve_qiime_environment(value, str)
                run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
