import importlib.util
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).parents[1] / "workflow" / "scripts" / "resolve_pr2_plastids.py"
)
SPEC = importlib.util.spec_from_file_location("resolve_pr2_plastids", SCRIPT)
RESOLVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RESOLVER)


def taxonomy_row(feature_id, taxonomy, confidence):
    return {
        "feature_id": feature_id,
        "taxonomy": taxonomy,
        "confidence_text": str(confidence),
        "confidence": confidence,
    }


class PR2ChloroplastResolutionTests(unittest.TestCase):
    def test_confident_pr2_plastids_override_silva_and_rescue_missed_asvs(self):
        silva = [
            taxonomy_row("confirmed", "d__Bacteria; p__Cyanobacteriota; o__Chloroplast", 0.99),
            taxonomy_row("rescued", "d__Bacteria; p__Minisyncoccota", 0.96),
            taxonomy_row("bacterium", "d__Bacteria; p__Pseudomonadota", 0.98),
        ]
        pr2 = [
            taxonomy_row("confirmed", "Eukaryota; Archaeplastida:plas; Chlorophyta", 0.95),
            taxonomy_row("rescued", "Eukaryota; Stramenopiles:plas; Ochrophyta", 0.91),
            taxonomy_row("bacterium", "Unassigned", 0.32),
        ]

        final_rows, audit_rows, summary = RESOLVER.resolve_taxonomies(
            silva, pr2, 0.7
        )
        final = {row["Feature ID"]: row["Taxon"] for row in final_rows}
        audit = {row["ASV_hash"]: row for row in audit_rows}

        self.assertIn(":plas", final["confirmed"])
        self.assertIn(":plas", final["rescued"])
        self.assertEqual(final["bacterium"], silva[2]["taxonomy"])
        self.assertEqual(
            audit["confirmed"]["Chloroplast_detection_source"],
            "PR2_confirmed_SILVA",
        )
        self.assertEqual(
            audit["rescued"]["Chloroplast_detection_source"], "PR2_rescue"
        )
        self.assertEqual(summary["PR2_rescued_ASVs"], 1)
        self.assertEqual(summary["PR2_confirmed_SILVA_ASVs"], 1)

    def test_low_confidence_pr2_plastid_does_not_override_silva(self):
        silva = [taxonomy_row("low", "d__Bacteria; p__Candidatus_Peribacteria", 0.88)]
        pr2 = [taxonomy_row("low", "Eukaryota; Alveolata:plas", 0.69)]

        final_rows, audit_rows, summary = RESOLVER.resolve_taxonomies(
            silva, pr2, 0.7
        )

        self.assertEqual(final_rows[0]["Taxon"], silva[0]["taxonomy"])
        self.assertEqual(
            audit_rows[0]["Chloroplast_detection_source"],
            "PR2_low_confidence_not_accepted",
        )
        self.assertEqual(summary["PR2_low_confidence_not_accepted_ASVs"], 1)

    def test_unconfirmed_silva_chloroplast_is_retained(self):
        silva = [
            taxonomy_row("silva-only", "d__Bacteria; p__Cyanobacteriota; o__Chloroplast", 0.93)
        ]
        pr2 = [taxonomy_row("silva-only", "Eukaryota; Opisthokonta", 0.91)]

        final_rows, audit_rows, summary = RESOLVER.resolve_taxonomies(
            silva, pr2, 0.7
        )

        self.assertEqual(final_rows[0]["Taxon"], silva[0]["taxonomy"])
        self.assertEqual(
            audit_rows[0]["Chloroplast_detection_source"],
            "SILVA_chloroplast_unconfirmed",
        )
        self.assertEqual(summary["SILVA_chloroplast_unconfirmed_ASVs"], 1)


if __name__ == "__main__":
    unittest.main()
