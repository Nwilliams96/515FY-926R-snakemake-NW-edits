rule prepare_bbsplit_db:
    input:
        file1=DATABASE_PREFIX + "bbsplit-db/SILVA_132_and_PR2_EUK.cdhit95pc.fasta", file2=DATABASE_PREFIX + "bbsplit-db/SILVA_132_PROK.cdhit95pc.fasta",
    output:
        directory(DATABASE_PREFIX + "bbsplit-db/EUK-PROK-bbsplit-db/")
    conda:
        "../envs/bbmap.yaml"
    log:
        "logs/bbsplit_db_prep.log"
    priority: 50
    shell:
        "bbsplit.sh build=1 ref={input.file1},{input.file2} path={output}"
