rule compute_target_features:
    """Compute the same feature set used for training, on the target assembly."""
    input:
        fasta="{assembly_dir}/{hap}_scaffolded.fa",
        female_fastqs=expand(
            "{raw_dir}/{acc}_1.fastq.gz",
            raw_dir=config["paths"]["raw_reads_dir"],
            acc=config["reads"]["female_comparison_accessions"],
        ),
        male_bam="{assembly_dir}/{hap}.hifi_depth.bam",
    output:
        features="{features_dir}/{hap}_features.csv",
    params:
        female_fastq_args=lambda wc, input: " ".join(f"--female-fastq {f}" for f in input.female_fastqs),
    shell:
        """
        python workflow/scripts/compute_target_features_cli.py \
            --fasta {input.fasta} \
            {params.female_fastq_args} \
            --male-bam {input.male_bam} \
            --out {output.features}
        """


rule repeatmask_candidates:
    input:
        fasta="{features_dir}/{hap}_y_candidates.fa",
    output:
        masked="{features_dir}/{hap}_y_candidates.fa.masked",
    threads: 8
    shell:
        "RepeatMasker -pa {threads} -species rodentia {input.fasta}"
