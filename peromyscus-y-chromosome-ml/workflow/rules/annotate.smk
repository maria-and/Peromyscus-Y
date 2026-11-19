rule align_rnaseq_evidence:
    input:
        fasta="{features_dir}/{hap}_y_candidates.fa.masked",
        rnaseq=expand(
            "{raw_dir}/{acc}_1.fastq.gz",
            raw_dir=config["paths"]["raw_reads_dir"],
            acc=config["rnaseq_evidence_accessions"],
        ),
    output:
        bam="{annotation_dir}/{hap}_rnaseq.bam",
    threads: 8
    shell:
        """
        hisat2-build {input.fasta} {wildcards.annotation_dir}/{wildcards.hap}_idx
        hisat2 -p {threads} -x {wildcards.annotation_dir}/{wildcards.hap}_idx \
            -U $(echo {input.rnaseq} | tr ' ' ',') | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule braker_annotate:
    input:
        fasta="{features_dir}/{hap}_y_candidates.fa.masked",
        bam="{annotation_dir}/{hap}_rnaseq.bam",
        proteins="config/known_y_proteins.fa",
    output:
        gff="{annotation_dir}/{hap}_braker/braker.gff3",
    threads: 16
    shell:
        """
        braker.pl --genome={input.fasta} --prot_seq={input.proteins} --bam={input.bam} \
            --softmasking --threads={threads} \
            --workingdir={wildcards.annotation_dir}/{wildcards.hap}_braker
        """
