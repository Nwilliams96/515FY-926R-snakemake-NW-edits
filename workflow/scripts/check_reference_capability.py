"""Fail before read processing if the selected RESCRIPt release is unsupported.

Only capability inspection is performed here. Acquisition and all scientific
reference preparation remain separate, visible Snakemake rules.

"""

import json
from pathlib import Path


def check_capability(family, version, output):
    import qiime2
    import rescript
    import sklearn
    from qiime2.sdk import PluginManager

    action_name = "get_silva_data" if family == "SILVA" else "get_pr2_data"
    action = PluginManager().plugins["rescript"].actions[action_name]
    supported_type = action.signature.parameters["version"].qiime_type
    if version not in supported_type:
        raise ValueError(
            f"Installed RESCRIPt {rescript.__version__} cannot download {family} {version}. "
            f"Supported versions: {supported_type}. Select a supported training release, "
            f"supply {family}_download_link for a compatible pre-trained classifier, "
            "or use an environment with support for the requested release."
        )
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(database=family, version=version,
                                   qiime2=qiime2.__version__, rescript=rescript.__version__,
                                   sklearn=sklearn.__version__), indent=2) + "\n")


if "snakemake" in globals():
    check_capability(snakemake.params.family, snakemake.params.version, snakemake.output[0])
