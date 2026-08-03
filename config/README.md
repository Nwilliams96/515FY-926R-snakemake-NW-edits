# Workflow configuration

Place the tutorial-generated configuration package in this directory. The
workflow loads `config/config.yml` and expects `config/samples.tsv` and
`config/bioanalyzer.tsv`.

## Database setup

Set `use_preexisting_databases: false` to let the workflow download and train
its databases under the repository-local `databases/` directory.

Set `use_preexisting_databases: true` to skip all database download and
training rules. In that mode, `database_dir` must point to a database root with
this structure:

```text
database_dir/
  bbsplit-db/
    EUK-PROK-bbsplit-db/
  classification/
    PR2/
      pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_<forward>_<reverse>_dereplicated_final_classifier_USE_ME.qza
    SILVA/
      silva-ssu-nr99-tax-dereplicated-sliced_<forward>_<reverse>_dereplicated_final_classifier_USE_ME.qza
```

The `<forward>` and `<reverse>` portions must match `fwdPrimer` and
`revPrimer` in `config.yml`.

## Internal standards

Set `use_internal_standards: false` to omit internal-standard correction.

When it is `true`, provide exactly three IDs under `intstds`, include matching
`<ID>_ng` columns in `samples.tsv`, and provide `config/internal_stds.tsv` with
these columns:

```text
internal_std_ID	rRNA_copy_number	genome_len_bp	full_16S_sequence
```

The workflow creates the required FASTA and BLAST database files from that
table automatically.

See `config.example.yml` for a complete configuration example.
