rule classify_y_candidates:
    input:
        model="{model_dir}/model.pkl",
        features="{features_dir}/{hap}_features.csv",
    output:
        classified="{features_dir}/{hap}_classified.csv",
    shell:
        """
        python workflow/scripts/classify_contigs.py \
            --model {input.model} \
            --features {input.features} \
            --out {output.classified}
        """


rule extract_y_candidate_fasta:
    input:
        classified="{features_dir}/{hap}_classified.csv",
        fasta="{assembly_dir}/{hap}_scaffolded.fa",
    output:
        y_fasta="{features_dir}/{hap}_y_candidates.fa",
    shell:
        """
        awk -F',' 'NR>1 && $NF=="True" {{print $1}}' {input.classified} > {output.y_fasta}.ids
        seqkit grep -f {output.y_fasta}.ids {input.fasta} > {output.y_fasta}
        """
