rule download_sra:
    """Download and extract a single SRA run accession.

    NOTE: fasterq-dump requires the actual RUN accession (SRR...), not an
    experiment accession (SRX...) — SRX resolution has to happen before this
    rule, e.g. via NCBI's efetch, not inside fasterq-dump itself.
    """
    output:
        r1="{raw_dir}/{accession}_1.fastq.gz",
    params:
        outdir=lambda wc: wc.raw_dir,
    threads: 4
    shell:
        """
        prefetch {wildcards.accession}
        fasterq-dump {wildcards.accession} --split-files -e {threads} -O {params.outdir}
        gzip -f {params.outdir}/{wildcards.accession}_*.fastq
        """
