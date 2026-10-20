#OSF_IDs = ["eux4r", "eahds"] #only for PROK / EUK splitting, other datbases available at repo
#NAMES = ["and_PR2_EUK", "PROK"]
#database_dir=config["database_dir"]

rule initialize_database_directories:
    output:
        marker=DATABASE_PREFIX + ".database_directories.ready"
    params:
        bbsplit=DATABASE_PREFIX + "bbsplit-db",
        silva=DATABASE_PREFIX + "classification/SILVA",
        pr2=DATABASE_PREFIX + "classification/PR2"
    priority: 51
    shell:
        "mkdir -p {params.bbsplit:q} {params.silva:q} {params.pr2:q} "
        "&& touch {output.marker:q}"

rule download_prok_db:
    input:
        rules.initialize_database_directories.output.marker
    output:
        temp(DATABASE_PREFIX + "bbsplit-db/SILVA_132_PROK.cdhit95pc.fasta"),
    log:
        "logs/aria2_bbsplit_prok_db_download.log",
    params:
        url="https://osf.io/eahds/download",
        # OSF's checksum pins the original reference. This downloader writes
        # the exact absolute path; aria2 --out incorrectly prefixed the CWD.
        checksum="46683398d8859aeef635fff49112c4d6453d4ebd70bf8b9fee006cb5b84273b1",
        checksum_algorithm="sha256",
    resources:
        mem_mb=1024,
        runtime=30,
    priority: 50
    script:
        "../scripts/download_verified_reference.py"

rule download_euk_db:
    input:
        rules.initialize_database_directories.output.marker
    output:
        temp(DATABASE_PREFIX + "bbsplit-db/SILVA_132_and_PR2_EUK.cdhit95pc.fasta"),
    log:
        "logs/aria2_bbsplit_euk_db_download.log",
    params:
        url="https://osf.io/eux4r/download",
        checksum="c4a396035271db516a0054f1d750d3081ca5a11831b44756e91590a79ed0f109",
        checksum_algorithm="sha256",
    resources:
        mem_mb=1024,
        runtime=30,
    priority: 50
    script:
        "../scripts/download_verified_reference.py"
