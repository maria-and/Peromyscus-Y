"""
Build a labeled training set from a reference genome with a known, resolved Y
chromosome (e.g. mouse GRCm39: chrY vs. autosomes).

Usage:
    python build_training_labels.py \
        --reference mouse_GRCm39.fa \
        --female-fastq female_reads.fastq.gz \
        --male-bam male_reads.sorted.bam \
        --out training_labels.csv

Labels: 1 = Y chromosome, 0 = autosome/X. Windows (not whole chromosomes) are
used as training examples so the classifier sees realistic within-chromosome
feature variance rather than one row per chromosome.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd
from Bio import SeqIO

from compute_features import (
    FemaleKmerSet,
    gc_content,
    female_absent_proportion,
    male_depth_from_bam,
    windows,
)

Y_CHROM_NAMES = {"Y", "chrY", "chrY_random"}


def build_labels(
    reference_fasta: str,
    female_fastq: str | None,
    male_bam: str | None,
    window_size: int,
    kmer_size: int,
    max_female_reads: int | None,
) -> pd.DataFrame:
    female_kmers = FemaleKmerSet(k=kmer_size)
    if female_fastq:
        print(f"Building female k-mer set from {female_fastq} ...", file=sys.stderr)
        female_kmers.add_fastq(female_fastq, max_reads=max_female_reads)
        print(f"  {len(female_kmers)} unique k-mers loaded", file=sys.stderr)

    rows = []
    for record in SeqIO.parse(reference_fasta, "fasta"):
        label = 1 if record.id in Y_CHROM_NAMES else 0
        seq = str(record.seq)
        for start, end in windows(len(seq), window_size):
            window_seq = seq[start:end]
            row = {
                "chrom": record.id,
                "start": start,
                "end": end,
                "length": end - start,
                "gc_content": gc_content(window_seq),
                "female_absent_prop": (
                    female_absent_proportion(window_seq, female_kmers)
                    if female_fastq
                    else None
                ),
                "male_depth": (
                    male_depth_from_bam(male_bam, record.id, start, end)
                    if male_bam
                    else None
                ),
                "label": label,
            }
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference", required=True, help="FASTA with a known Y chromosome"
    )
    parser.add_argument(
        "--female-fastq", default=None, help="Female short-read FASTQ(.gz)"
    )
    parser.add_argument(
        "--male-bam", default=None, help="Male reads aligned to --reference"
    )
    parser.add_argument("--window-size", type=int, default=10_000)
    parser.add_argument("--kmer-size", type=int, default=25)
    parser.add_argument("--max-female-reads", type=int, default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df = build_labels(
        args.reference,
        args.female_fastq,
        args.male_bam,
        args.window_size,
        args.kmer_size,
        args.max_female_reads,
    )
    df.to_csv(args.out, index=False)
    print(
        f"Wrote {len(df)} labeled windows ({df['label'].sum()} Y, "
        f"{(df['label'] == 0).sum()} non-Y) to {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
