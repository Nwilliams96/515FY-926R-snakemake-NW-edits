rule cutadapt:
    input:
        r1=lambda wildcards: os.path.join(
            config["rawdatadir"], wildcards.sample + config["R1file_ending"]
        ),
        r2=lambda wildcards: os.path.join(
            config["rawdatadir"], wildcards.sample + config["R2file_ending"]
        ),
        # Check native database support before processing FASTQs. The pinned
        # wrapper requires exactly two total inputs, so its unchanged command
        # is explicit below and passes only r1/r2 to Cutadapt.
        database_checks=DATABASE_SETUP_CHECKS,
    output:
        fastq1=temp("results/00-trimmed/{sample}.1.fastq"),
        fastq2=temp("results/00-trimmed/{sample}.2.fastq"),
        qc="results/00-trimmed/{sample}.qc.txt",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-g ^" + config["fwdPrimer"] + " -G ^" + config["revPrimer"],
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--no-indels --pair-filter=any --error-rate=0.2 --discard-untrimmed",
    log:
        "logs/00-trimming/{sample}.log",
    threads: 4  # set desired number of threads here
    conda:
        "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v7.0.0/bio/cutadapt/pe/environment.yaml"
    shell:
        "cutadapt --cores {threads} {params.adapters} {params.extra} "
        "-o {output.fastq1:q} -p {output.fastq2:q} "
        "{input.r1:q} {input.r2:q} > {output.qc:q} 2> {log:q}"

rule compress_trimmed_output_r1:
    input:
        "results/00-trimmed/{sample}.1.fastq",
    output:
        "results/00-trimmed/{sample}.1.fastq.gz",
    threads: 1
    log:
        "logs/bgzip/{sample}_compress-trimmed-r1.log",
    wrapper:
        "v7.2.0/bio/bgzip"

rule compress_trimmed_output_r2:
    input:
        "results/00-trimmed/{sample}.2.fastq"
    output:
        "results/00-trimmed/{sample}.2.fastq.gz"
    threads: 1
    log:
        "logs/bgzip/{sample}_compress-trimmed-r2.log",
    wrapper:
        "v7.2.0/bio/bgzip"
