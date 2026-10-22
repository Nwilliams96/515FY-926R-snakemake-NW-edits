"""Check the website's ordered-list ISD configuration and downstream names."""
from itertools import combinations
from pathlib import Path
import re
import unittest

from snakemake.exceptions import WorkflowError

ROOT = Path(__file__).resolve().parents[1]
COMMON = (ROOT / "workflow/rules/common.smk").read_text()


class InternalStandardConfigurationTests(unittest.TestCase):
    def resolve(self, config):
        block = COMMON[COMMON.index("configured_intstds ="):COMMON.index("NEW_AMPLICON_CONCENTRATIONS_FILE =")]
        scope = dict(config=dict(studyName="fixture", **config),
                     USE_INTERNAL_STANDARDS=config.get("use_internal_standards", False),
                     WorkflowError=WorkflowError, re=re, combinations=combinations)
        exec(compile(block, "common.smk:internal-standards", "exec"), scope)
        return scope

    def test_current_list_configs_preserve_ids_order_and_named_outputs(self):
        for ids in (["BP"], ["BP", "DR"], ["BP", "DR", "TT"],
                    ["Custom-A", "Custom_B", "SpikeIn3", "SpikeIn4"]):
            with self.subTest(ids=ids):
                result = self.resolve(dict(use_internal_standards=True, intstds=ids))
                self.assertEqual(result["INTERNAL_STANDARD_IDS"], ids)
                self.assertEqual(result["INTERNAL_STANDARD_SLOTS"],
                                 [f"isd_{i}" for i in range(1, len(ids) + 1)])
                self.assertEqual(result["INTERNAL_STANDARD_ID_PAIRS"], list(combinations(ids, 2)))
                self.assertEqual(result["ALL_INTERNAL_STANDARD_IDS_STEM"], "_".join(ids))
                self.assertEqual(result["ISD_CORRECTED_LONG_TABLE"],
                                 "results/05-internal-std-corrected/fixture."
                                 + "_".join(ids) + ".ISD_corrected_long_data.tsv")
                self.assertEqual(result["INTERNAL_STANDARD_METHOD_STEMS"],
                                 [f"{id}_recovery_ratio" for id in ids]
                                 + [f"mean_{'_'.join(ids)}_recovery_ratio",
                                    f"median_{'_'.join(ids)}_recovery_ratio"]
                                 + [f"mean_{a}_and_{b}_recovery_ratio" for a,b in combinations(ids,2)])

    def test_list_entries_are_trimmed_and_converted_to_text(self):
        result = self.resolve(dict(use_internal_standards=True, intstds=[" BP ", "\tDR", 123]))
        self.assertEqual(result["INTERNAL_STANDARD_IDS"], ["BP", "DR", "123"])

    def test_missing_or_empty_list_remains_valid_when_standards_are_disabled(self):
        for config in ({}, {"use_internal_standards": False},
                       {"use_internal_standards": False, "intstds": []}):
            with self.subTest(config=config):
                self.assertEqual(self.resolve(config)["INTERNAL_STANDARD_IDS"], [])

    def test_legacy_mappings_and_other_nonlist_types_fail_clearly(self):
        for value in ({"intstd1": "BP", "intstd2": "DR"}, {}, "BP", None, 1, ("BP",)):
            for enabled in (True, False):
                with self.subTest(value=value, enabled=enabled), self.assertRaisesRegex(
                    WorkflowError, "intstds must be an ordered list of internal-standard IDs"
                ):
                    self.resolve(dict(use_internal_standards=enabled, intstds=value))

    def test_enabled_standards_still_require_nonempty_unique_safe_ids(self):
        for config, message in (
            ({"use_internal_standards": True}, "does not contain any IDs"),
            ({"use_internal_standards": True, "intstds": []}, "does not contain any IDs"),
            ({"use_internal_standards": True, "intstds": [""]}, "may contain only"),
            ({"use_internal_standards": True, "intstds": ["bad name"]}, "may contain only"),
            ({"use_internal_standards": True, "intstds": ["BP", " BP "]}, "must be unique"),
        ):
            with self.subTest(config=config), self.assertRaisesRegex(WorkflowError, message):
                self.resolve(config)


if __name__ == "__main__":
    unittest.main()
