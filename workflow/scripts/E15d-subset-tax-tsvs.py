"""Create metazoan and non-metazoan TSVs without requiring non-empty QIIME tables."""

from pathlib import Path


SILVA_METAZOA = (
    "p__Annelida",
    "p__Apicomplexa",
    "p__Arthropoda",
    "p__Cnidaria",
    "p__Ctenophora",
    "p__Echinodermata",
    "p__Holozoa",
    "p__Mollusca",
    "p__Porifera",
    "p__Tunicata",
    "p__Vertebrata",
)


def split_table(source, included, excluded, taxonomy_matches):
    lines = Path(source).read_text().splitlines(keepends=True)
    if len(lines) < 2 or not lines[1].startswith("#OTU ID\t"):
        raise ValueError(f"Unexpected BIOM TSV header in {source}")

    header = lines[:2]
    included_rows = []
    excluded_rows = []
    for line in lines[2:]:
        taxonomy = line.rstrip("\r\n").split("\t")[-1]
        (included_rows if taxonomy_matches(taxonomy) else excluded_rows).append(line)

    Path(included).write_text("".join(header + included_rows))
    Path(excluded).write_text("".join(header + excluded_rows))


split_table(
    snakemake.input.all18SSILVA,
    snakemake.output.includemetazoaSILVAtablebiomtaxtsv,
    snakemake.output.excludemetazoaSILVAtablebiomtaxtsv,
    lambda taxonomy: any(term in taxonomy for term in SILVA_METAZOA),
)

split_table(
    snakemake.input.all18SPR2,
    snakemake.output.includemetazoaPR2tablebiomtaxtsv,
    snakemake.output.excludemetazoaPR2tablebiomtaxtsv,
    lambda taxonomy: "Metazoa" in taxonomy,
)
