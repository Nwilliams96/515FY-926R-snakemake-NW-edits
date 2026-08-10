rule prepare_bbsplit_db:
    input:
        file1=DATABASE_PREFIX + "bbsplit-db/SILVA_132_and_PR2_EUK.cdhit95pc.fasta", file2=DATABASE_PREFIX + "bbsplit-db/SILVA_132_PROK.cdhit95pc.fasta",
    output:
        directory(DATABASE_PREFIX + "bbsplit-db/EUK-PROK-bbsplit-db/")
    conda:
        "../envs/bbmap.yaml"
    log:
        "logs/bbsplit_db_prep.log"
    threads: 8
    resources:
        mem_mb=64000,
        runtime=240,
    priority: 50
    shell:
        "bbsplit.sh build=1 threads={threads} ref={input.file1:q},{input.file2:q} "
        "path={output:q} 2> {log:q}"
