# Code provenance and human-review ledger

This document records significant AI-assisted changes to the executable 515Y–926R eASV workflow. It is intended to make scientific and engineering review easier; it is not a substitute for reviewing the code, tests, and representative outputs.

Git records these commits under Nathan Williams's configured Git identity. Git authorship alone does not distinguish text written directly by Nathan from code drafted, revised, or debugged with Codex. The entries below are therefore reconstructed from the project conversation history and repository history. For a conservative audit, executable workflow changes made after the initial import ([`1644b19`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/1644b19)) are treated as AI-assisted unless specifically documented otherwise.

All linked historical commits have corresponding objects in the local Git history. The link for the new local refactor commit `1a25e17` will become available on GitHub only if this branch is later pushed; this implementation intentionally does not push it without separate authorization.

Review values are `Pending`, `Accepted`, `Revise`, or `Reverted`. All major entries begin as `Pending`; a named human reviewer should update the record only after examining the diff and relevant outputs.

## Scientific and workflow invariants

- Canonical prokaryotic ASVs are classified with SILVA 144, then screened independently with PR2 5.1.1 for plastid recovery. A PR2 assignment replaces SILVA only when the lineage contains `:plas` and confidence is at least 0.7.
- The official SILVA 144 uniform classifier stops at genus. Species labels are provisional exact matches from the matching SILVA DADA2 reference and are accepted only when the exact-match genus agrees with the classifier genus.
- Canonical merged 16S features remain separate from the linked/unmerged 16S sensitivity supplement.
- Eukaryotic reads are fused into the established concatenated 18S representation and classified once with PR2.
- Genomic internal standards may be 16S or 18S and may number one or more. Detected standard ASVs are excluded from ecological tables and report composition plots.
- Shared databases are release-pinned and may be reused, but incomplete files must never be published as completed shared resources.

## Major design and scientific changes

### [`18ead1d`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/18ead1d) — Optional database setup

- **Files:** `workflow/Snakefile`, database rules, configuration handling.
- **Request/problem:** Avoid downloading and training databases when compatible resources already exist.
- **Nathan's rationale:** Shared laboratory resources should be reusable across study-specific clones.
- **Codex rationale:** Separate configuration-controlled database branches from sample processing while preserving explicit dependencies.
- **Behavioral impact:** Database preparation can be skipped in favor of configured resources.
- **Tradeoffs:** Incorrectly marked pre-existing resources now fail validation or appear as missing inputs rather than being silently rebuilt.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`054cb57`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/054cb57) — Databases outside run clones

- **Files:** shared database paths and related configuration.
- **Request/problem:** A fresh clone is created for each project, so large databases must not be duplicated inside it.
- **Nathan's rationale:** Keep projects organized without repeatedly storing the same references.
- **Codex rationale:** Resolve all database products beneath one configured external root.
- **Behavioral impact:** Database files are shared by project clones.
- **Tradeoffs:** Concurrent workflows now require locking and atomic publication, added later.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`069addf`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/069addf) — Self-contained pipeline report

- **Files:** report rule and `workflow/scripts/generate_pipeline_report.py`.
- **Request/problem:** Provide a concise HTML account of configuration, QC, read retention, composition, taxonomy, and internal-standard results.
- **Nathan's rationale:** Users need a rapid first inspection before detailed analysis.
- **Codex rationale:** Generate a portable, dependency-free HTML report from tracked pipeline outputs.
- **Behavioral impact:** A report becomes a final workflow target and export.
- **Tradeoffs:** Report aggregation adds code and must not be mistaken for inferential analysis.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`3d56c58`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/3d56c58) — Configurable internal-standard identities

- **Files:** internal-standard rules, R correction code, output naming.
- **Request/problem:** Replace hard-coded BP/DR/TT labels with configurable standard identities.
- **Nathan's rationale:** Other users employ different genomic standards and need recognizable figure/output labels.
- **Codex rationale:** Pass stable standard IDs from configuration through detection, correction, plots, and filenames.
- **Behavioral impact:** Standard identities are no longer laboratory-specific constants.
- **Tradeoffs:** IDs must be unique and safe for paths/column names.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`bdc17ea`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/bdc17ea) — Variable internal-standard count

- **Files:** `workflow/rules/common.smk`, internal-standard rules, correction script.
- **Request/problem:** Support datasets with one, two, three, or more standards.
- **Nathan's rationale:** The number of standards varies by experiment.
- **Codex rationale:** Derive standard slots, pairwise methods, required columns, and outputs from an ordered list.
- **Behavioral impact:** The optional correction branch scales dynamically with the configured list.
- **Tradeoffs:** Combination outputs grow with the number of standards and require explicit naming.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`f7aba4f`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/f7aba4f) — Configurable DADA2 parameters

- **Files:** prokaryotic/eukaryotic DADA2 rules and shell scripts, configuration/report plumbing.
- **Request/problem:** Permit controlled tests of expected-error, truncation, merging, pooling, chimera, and learning settings.
- **Nathan's rationale:** Evaluate read/ASV recovery across datasets without editing scripts.
- **Codex rationale:** Preserve historical defaults while recording the effective parameters in configuration and reports.
- **Behavioral impact:** DADA2 behavior is study-configurable and auditable.
- **Tradeoffs:** Parameter changes can make comparisons invalid unless baselines and changed factors are recorded.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`86ca36d`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/86ca36d) — SILVA 144 adoption

- **Files:** database constants/rules, QIIME environment, taxonomy parsing.
- **Request/problem:** Make SILVA 144 the default prokaryotic reference.
- **Nathan's rationale:** Use the current release and its 515F/926R classifier.
- **Codex rationale:** Pin URL, checksum, classifier/environment compatibility, and rank-aware parsing so the new Kingdom rank cannot shift downstream columns.
- **Behavioral impact:** Canonical prokaryotic taxonomy uses SILVA 144 and QIIME 2 2026.7.
- **Tradeoffs:** Results are not directly identical to older SILVA releases; classifier filenames and environment compatibility are release-specific.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`c7a1009`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/c7a1009) — Provisional exact species matching

- **Files:** `workflow/scripts/assign_silva_species.R` and its rule/output summary.
- **Request/problem:** Preserve species information not emitted by the uniform genus-level SILVA classifier.
- **Nathan's rationale:** Species-level labels are useful when the reference contains an exact sequence match.
- **Codex rationale:** Accept only unambiguous exact matches whose genus agrees with the classifier; otherwise leave species blank.
- **Behavioral impact:** Some SILVA lineages receive an appended `s__Genus species` label and an audit summary.
- **Tradeoffs:** SILVA organism names are not curated species calls; matches remain provisional and can be memory intensive, hence later batching.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`92f67ee`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/92f67ee) — PR2 plastid resolution

- **Files:** prokaryotic classification/subsetting rules, reconciliation script, category split logic.
- **Request/problem:** SILVA 144 assigned dominant plastid-derived 16S ASVs as cyanobacterial lineages, losing useful chloroplast identities.
- **Nathan's rationale:** Recover chloroplast sequences using PR2 before chloroplast/cyanobacteria outputs are split.
- **Codex rationale:** Classify the same canonical ASVs independently, retain SILVA by default, and accept only explicit confident PR2 `:plas` calls with an audit trail.
- **Behavioral impact:** Confirmed and rescued plastids receive PR2 lineages before category splitting.
- **Tradeoffs:** Adds a second classifier pass and an explicit 0.7 decision threshold requiring scientific review.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`f13649a`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/f13649a) — Linked/unmerged 16S supplement

- **Files:** second DADA2 path, classification/export rules, supplementary-table script.
- **Request/problem:** Examine paired reads lost because canonical 16S pairs do not merge at the configured truncation lengths.
- **Nathan's rationale:** Retain potentially informative features for sensitivity analysis without altering the main table.
- **Codex rationale:** Run DADA2 with linked-pair retention, subtract features already present in the canonical run, classify the remainder, and export separately.
- **Behavioral impact:** Two supplementary TSVs are produced automatically and copied to Results-Export.
- **Tradeoffs:** Linked pairs are not directly interchangeable with merged ASVs and must not be automatically combined or used for canonical richness/abundance.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`8d63a81`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/8d63a81) — Shared PR2 concurrency safety

- **Files:** PR2 construction workflow and shared-publication logic.
- **Request/problem:** Simultaneous study runs deleted or replaced shared PR2 intermediates used by other jobs.
- **Nathan's rationale:** Multiple datasets should run concurrently against one shared database directory.
- **Codex rationale:** Build privately, validate, serialize publication with a filesystem lock, and atomically replace only the completed classifier.
- **Behavioral impact:** Partial PR2 artifacts are no longer exposed as shared completed outputs.
- **Tradeoffs:** The first implementation consolidated scientific steps in Python, reducing DAG visibility; that design is reversed by `1a25e17` while retaining the safety properties.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`f380a1b`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/f380a1b) — Concatenated 18S internal standards

- **Files:** internal-standard FASTA preparation/detection and tests.
- **Request/problem:** Whole-genome 18S standards need to match the pipeline's R1-plus-reverse-complemented-R2 concatenated representation.
- **Nathan's rationale:** Permit genomic internal standards from either marker branch without requiring users to pre-transform sequences.
- **Codex rationale:** Derive the expected amplicon representation from the full SSU sequence and configured primers/truncation, then search both branches.
- **Behavioral impact:** 16S and 18S genomic standards can be detected with the same configuration model.
- **Tradeoffs:** Primer orientation and complete reference sequence quality are now critical inputs.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

### [`1a25e17`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/1a25e17) — Transparent classification rules

- **Files:** `workflow/rules/prepare_classification_dbs.smk`, `workflow/rules/02-denoise-and-export-prok.smk`, small Python transformations/publication helper, tests.
- **Request/problem:** Consolidated Python obscured PR2 construction and QIIME/BIOM export operations from the Snakemake DAG.
- **Nathan's rationale:** Restore the original teaching-friendly, rule-per-step overview so collaborators can understand and modify each operation.
- **Codex rationale:** Put scientific and file-conversion commands back into named rules while limiting Python to pure TSV transformations and locked atomic publication.
- **Behavioral impact:** PR2 download through classifier training, chloroplast reconciliation, and linked-supplement exports are visible in the DAG; filenames, thresholds, and final outputs remain unchanged.
- **Tradeoffs:** The rule file is longer, and independently launched projects may duplicate private PR2 construction work before one publishes; shared outputs remain protected.
- **Status:** Pending
- **Reviewer:**
- **Review date:**
- **Notes:**

## Operational fixes

These changes address execution, compatibility, metadata, memory, or empty-output behavior rather than intentionally changing the scientific interpretation. Each link opens the exact GitHub diff.

| Commit | Operational purpose |
| --- | --- |
| [`3e883d4`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/3e883d4) | Point the CARC launcher at the workflow Snakefile. |
| [`691e2a0`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/691e2a0) | Require the generated study configuration. |
| [`90ea921`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/90ea921) | Keep tutorial support files inside `config/`. |
| [`caac19b`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/caac19b) | Replace the original CARC-specific generic launcher. |
| [`8f4dd0b`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/8f4dd0b) | Clarify database inputs and remove invalid directory flags. |
| [`719c12d`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/719c12d) | Add the explicitly USC-only CARC launcher. |
| [`56f0941`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/56f0941) | Fix amplicon-molarity input handling and empty eukaryotic subsets. |
| [`1d8b076`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/1d8b076) | Merge the completed development branch into the then-current workflow. |
| [`c9015bc`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/c9015bc) | Preserve textual sample-type values. |
| [`712861c`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/712861c) | Remove a remote BBTools wrapper dependency. |
| [`a99d5d5`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/a99d5d5) | Align CARC CPU/memory allocation with rule scheduling. |
| [`9d08f79`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/9d08f79), [`1f8aec8`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/1f8aec8) | Label CARC jobs with the study name. |
| [`dad64a7`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/dad64a7) | Repair future-dated transferred config timestamps. |
| [`0410a1f`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/0410a1f) | Reuse completed Conda rule environments across clones. |
| [`7add214`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/7add214) | Isolate each Slurm job's Conda package cache to avoid stale NFS handles. |
| [`72c8c7c`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/72c8c7c) | Configure BBsplit, SILVA, and PR2 pre-existing-resource choices independently. |
| [`a203833`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/a203833) | Handle empty ProPortal classifications. |
| [`59b407d`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/59b407d) | Cap concurrent BBTools Java memory. |
| [`dafd00d`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/dafd00d) | Read numeric sample identifiers as text. |
| [`07fead1`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/07fead1) | Accept SILVA 144 cyanobacterial naming in subsets. |
| [`7ab256f`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/7ab256f) | Remove unused Deblur dependencies from the QIIME environment. |
| [`29d757d`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/29d757d) | Correct PR2 download placement. |
| [`53de9ae`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/53de9ae) | Support QIIME 2 PR2 extraction-statistics output. |
| [`5c6c1aa`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/5c6c1aa) | Support QIIME 2 DADA2 transition-statistics output. |
| [`3dc9610`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/3dc9610) | Exclude taxonomy metadata from the count pivot. |
| [`bae067b`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/bae067b), [`cd392d2`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/cd392d2), [`46ef33c`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/46ef33c), [`0dd4e51`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/0dd4e51) | Record test-branch truncation changes and restore historical defaults on the relevant branches. |
| [`8631e44`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/8631e44) | Merge the chloroplast and linked-read test work. |
| [`d1f7d7d`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/d1f7d7d) | Make the separate linked-read supplement an automatic output. |
| [`9e3a185`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/9e3a185) | Bound PR2 and exact-species memory use. |
| [`d230025`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/d230025) | Lock, checksum, and atomically publish shared SILVA references. |
| [`55ea31a`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/55ea31a) | Remove duplicate 18S PR2 classification and relabel outputs honestly. |

## Reporting and user-interface changes

The following commits primarily change report presentation, summaries, or user-facing packaging rather than the underlying denoising/classification computation:

- [`d97329b`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/d97329b) fixed report execution under Snakemake.
- [`b73f253`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/b73f253) collected final deliverables in Results-Export.
- [`090cef2`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/090cef2) added stage-specific DADA2 loss summaries.
- [`895ca92`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/895ca92) removed confusing quality/loss sections.
- [`5f3378e`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/5f3378e) added sequence-assignment summaries.
- [`75cb2d1`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/75cb2d1) reported effective DADA2 parameters.
- [`51d7855`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/51d7855) added interactive taxonomy metadata filters.
- [`2dc44ab`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/2dc44ab) expanded report inputs and pipeline configuration coverage.
- [`4541685`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/4541685) added the QIIME-style taxonomy bar plot.
- [`8fb177b`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/8fb177b) improved ISD exports, reports, and packaging.
- [`a1eef28`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/a1eef28) moved taxonomy labels below bars.
- [`78d046a`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/78d046a) handled missing values in charts.
- [`52677ec`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/52677ec) improved axes, palettes, and read-loss displays.
- [`da22702`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/da22702) removed detected ISD ASVs from final biological tables and report plots.
- [`77c91b8`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/77c91b8) displayed exact 16S/18S correction factors.
- [`92c6c58`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/92c6c58) embedded QIIME quality profiles.

## Documentation-only and repository-administration changes

These commits did not intentionally revise executable scientific behavior:

- [`a953109`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/a953109), [`74d56f0`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/74d56f0), and [`aa76c02`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/aa76c02) updated requirements and README material.
- [`383a8dc`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/383a8dc) linked the public tutorial.
- [`0f11aa1`](https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/commit/0f11aa1) corrected the documented CARC submission command.

## Why the consolidated PR2 builder was replaced

Commit `8d63a81` correctly addressed a real concurrency failure: independent Snakemake processes were modifying the same shared intermediate PR2 files. Its implementation moved the full build into `ensure_pr2_classifier.py`, where locking and private staging were straightforward but the scientific steps disappeared from the Snakemake DAG.

Commit `1a25e17` restores visible rules for download, decompression, FASTA/taxonomy reformatting, imports, culling, full-reference dereplication, primer extraction, extracted-reference dereplication, classifier training, validation, and publication. Intermediate artifacts now live in each project's `.snakemake/database-build/PR2/` directory. The remaining Python publication helper performs only QIIME validation, lock acquisition, and atomic publication of a completed classifier. Thus the safety fix is retained without hiding the scientific workflow.

## BBsplit teaching-run download incident — unresolved

A teaching run reportedly timed out in `download_prok_db`, and a possible failure to create the BBsplit directory was suggested. No downloader change is made here because the required failure evidence has not been provided. The current `initialize_database_directories` rule creates the BBsplit, SILVA, and PR2 directories in the same command; available evidence therefore does not support a selective BBsplit directory-creation failure.

Before a separate approved patch, review:

- `logs/aria2_bbsplit_prok_db_download.log`
- the corresponding Snakemake or Slurm error log
- destination permissions and filesystem state
- whether the failure was an OSF timeout, redirect/wrapper retrieval problem, retry policy, permission problem, or output-path mismatch

## Review checklist

- [ ] Inspect each major linked diff and update its review fields.
- [ ] Confirm SILVA 144, PR2 5.1.1, QIIME 2 2026.7, and the 0.7 plastid threshold are scientifically intended.
- [ ] Review exact-species limitations and wording before publication.
- [ ] Compare resolved taxonomy, routing audit/summary, and linked-supplement fixtures across the transparent refactor.
- [ ] Dry-run both missing-PR2 and valid-pre-existing-PR2 configurations in an environment with Snakemake.
- [ ] Run a representative CARC dataset only after explicit authorization.
- [ ] Resolve the BBsplit incident only after the requested logs are available.

