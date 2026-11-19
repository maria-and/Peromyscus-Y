rule hifiasm_assembly:
    """Haplotype-resolved assembly from HiFi + Hi-C/Omni-C reads."""
    input:
        hifi="{raw_dir}/all_hifi.fastq.gz",
        hic_r1="{raw_dir}/hic_1.fastq.gz",
        hic_r2="{raw_dir}/hic_2.fastq.gz",
    output:
        hap1="{assembly_dir}/asm.hic.hap1.p_ctg.gfa",
        hap2="{assembly_dir}/asm.hic.hap2.p_ctg.gfa",
    threads: config["assembly"]["hifiasm_threads"]
    resources:
        mem_mb=config["assembly"]["hifiasm_mem_gb"] * 1024,
    shell:
        """
        hifiasm -o {wildcards.assembly_dir}/asm -t {threads} \
            --h1 {input.hic_r1} --h2 {input.hic_r2} {input.hifi}
        """


rule gfa_to_fasta:
    input:
        gfa="{assembly_dir}/asm.hic.{hap}.p_ctg.gfa",
    output:
        fasta="{assembly_dir}/{hap}.fa",
    shell:
        """awk '/^S/{{print ">"$2"\\n"$3}}' {input.gfa} > {output.fasta}"""


rule scaffold_with_hic:
    """Scaffold a haplotype assembly using Hi-C/Omni-C contact data (YaHS)."""
    input:
        fasta="{assembly_dir}/{hap}.fa",
        hic_bam="{assembly_dir}/{hap}.hic.bam",
    output:
        scaffolded="{assembly_dir}/{hap}_scaffolded.fa",
    shell:
        "yahs {input.fasta} {input.hic_bam} -o {wildcards.assembly_dir}/{wildcards.hap}_scaffolded"
