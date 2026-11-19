"""
CLI wrapper around compute_contig_features() for use in the Snakemake pipeline
(kept as a standalone script, like the rest of workflow/scripts, rather than a
Snakemake `script:` block, so it can also be run/tested outside Snakemake).

Usage:
    python compute_target_features_cli.py \
        --fasta hap1_scaffolded.fa \
        --female-fastq female1.fastq.gz --female-fastq female2.fastq.gz \
        --male-bam hap1.hifi_depth.bam \
        --repeatmasker-out hap1.fa.out \
        --kmer-size 25 \
        --out hap1_features.csv
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from compute_features import FemaleKmerSet, compute_contig_features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument(
        "--female-fastq", action="append", default=[], help="repeatable"
    )
    parser.add_argument("--male-bam", default=None)
    parser.add_argument("--repeatmasker-out", default=None)
    parser.add_argument("--kmer-size", type=int, default=25)
    parser.add_argument("--max-female-reads-per-file", type=int, default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    female_kmers = FemaleKmerSet(k=args.kmer_size)
    for fq in args.female_fastq:
        print(f"Loading female k-mers from {fq} ...", file=sys.stderr)
        female_kmers.add_fastq(fq, max_reads=args.max_female_reads_per_file)
    print(f"Female k-mer set size: {len(female_kmers)}", file=sys.stderr)

    rows = compute_contig_features(
        args.fasta,
        female_kmers,
        bam_path=args.male_bam,
        repeatmasker_out=args.repeatmasker_out,
    )
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"Wrote features for {len(rows)} contigs to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
