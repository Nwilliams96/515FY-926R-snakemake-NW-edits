"""Copy the small, user-facing deliverables into Results-Export."""

from pathlib import Path
import shutil


def export_results(sources, destinations, readme_path, study_name):
    if len(sources) != len(destinations):
        raise ValueError("Every export source must have one destination")

    copied_names = []
    for source, destination in zip(sources, destinations):
        source_path = Path(source)
        destination_path = Path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination_path)
        copied_names.append(destination_path.name)

    readme = Path(readme_path)
    readme.parent.mkdir(parents=True, exist_ok=True)
    contents = [
        f"Results export for {study_name}",
        "",
        "This folder contains the small, user-facing outputs from the completed pipeline run.",
        "The HTML report is self-contained and can be opened directly in a web browser.",
        "",
        "Files:",
        *[f"- {name}" for name in copied_names],
        "- README.txt",
        "",
        "Keep the complete results/ directory if detailed QIIME 2 artifacts or intermediate files are needed.",
        "",
    ]
    readme.write_text("\n".join(contents), encoding="utf-8")


if "snakemake" in globals():
    export_results(
        list(snakemake.input.data_files),
        list(snakemake.output.data_files),
        str(snakemake.output.readme),
        str(snakemake.config["studyName"]),
    )
