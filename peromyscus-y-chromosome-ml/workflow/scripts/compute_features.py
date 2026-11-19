"""
Per-contig / per-window feature extraction for Y-candidate classification.

Features computed:
  - gc_content: fraction of G/C bases
  - female_absent_prop: proportion of this sequence's k-mers that are ABSENT
    from a female k-mer set (high = Y-like, since females have no Y)
  - male_depth: median male read depth across the window/contig
  - repeat_frac: fraction of bases masked by RepeatMasker (if a .out file is given)

Designed to be called both when building labeled training data (from a reference
with a known Y) and when scoring a real target assembly's contigs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

# pysam and biopython are only needed for the BAM/FASTA-parsing helpers below —
# imported lazily inside those functions so the pure feature-math (gc_content,
# kmerize, female_absent_proportion, windows) stays importable/testable without
# requiring those (heavier, C-extension-backed) dependencies installed.


def gc_content(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    denom = gc + at
    return gc / denom if denom else 0.0


def kmerize(seq: str, k: int) -> set[str]:
    seq = seq.upper()
    if len(seq) < k:
        return set()
    return {
        seq[i : i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i : i + k]
    }


class FemaleKmerSet:
    """Lightweight k-mer presence/absence lookup built from female reads.

    Backed by a Python set for clarity/testability at small scale; for a real
    whole-genome female read set, swap this for a Bloom filter (e.g.
    `pybloomfiltermmap3`) or a disk-backed k-mer counter (DSK, KMC) — the
    interface (`__contains__`) stays the same either way.
    """

    def __init__(self, k: int = 25):
        self.k = k
        self._kmers: set[str] = set()

    def add_sequence(self, seq: str) -> None:
        self._kmers.update(kmerize(seq, self.k))

    def add_fastq(self, fastq_path: str | Path, max_reads: int | None = None) -> None:
        import gzip

        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        with opener(fastq_path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:  # sequence line of a FASTQ record
                    self.add_sequence(line.strip())
                if max_reads is not None and i // 4 >= max_reads:
                    break

    def __contains__(self, kmer: str) -> bool:
        return kmer in self._kmers

    def __len__(self) -> int:
        return len(self._kmers)


def female_absent_proportion(seq: str, female_kmers: FemaleKmerSet) -> float:
    """Proportion of this sequence's k-mers that do NOT appear in the female set.

    High values (close to 1.0) mean this sequence looks Y-like: almost none of
    its k-mers exist in a female genome's reads at all. Low values mean the
    sequence is mostly shared with the female genome (i.e. autosomal/X-like).
    """
    kmers = kmerize(seq, female_kmers.k)
    if not kmers:
        return 0.0
    absent = sum(1 for km in kmers if km not in female_kmers)
    return absent / len(kmers)


def male_depth_from_bam(
    bam_path: str | Path, contig: str, start: int, end: int
) -> float:
    """Median male read depth over [start, end) on `contig`, via pysam."""
    import numpy as np
    import pysam

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        depths = [
            col.nsegments
            for col in bam.pileup(contig, start, end, truncate=True, min_base_quality=0)
        ]
    return float(np.median(depths)) if depths else 0.0


def repeat_fraction_from_repeatmasker(
    out_path: str | Path, contig: str, start: int, end: int
) -> float:
    """Fraction of [start, end) on `contig` masked, parsed from a RepeatMasker .out file."""
    masked = 0
    window_len = end - start
    with open(out_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 7 or not parts[0].isdigit():
                continue
            rm_contig, rm_start, rm_end = parts[4], int(parts[5]), int(parts[6])
            if rm_contig != contig:
                continue
            overlap_start = max(start, rm_start)
            overlap_end = min(end, rm_end)
            if overlap_end > overlap_start:
                masked += overlap_end - overlap_start
    return masked / window_len if window_len else 0.0


def windows(contig_length: int, window_size: int) -> Iterable[tuple[int, int]]:
    for start in range(0, contig_length, window_size):
        yield start, min(start + window_size, contig_length)


def compute_contig_features(
    fasta_path: str | Path,
    female_kmers: FemaleKmerSet,
    bam_path: str | Path | None = None,
    repeatmasker_out: str | Path | None = None,
) -> list[dict]:
    """Compute whole-contig-level features for every sequence in a FASTA file."""
    from Bio import SeqIO

    rows = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(record.seq)
        row = {
            "contig": record.id,
            "length": len(seq),
            "gc_content": gc_content(seq),
            "female_absent_prop": female_absent_proportion(seq, female_kmers),
        }
        if bam_path is not None:
            row["male_depth"] = male_depth_from_bam(bam_path, record.id, 0, len(seq))
        if repeatmasker_out is not None:
            row["repeat_frac"] = repeat_fraction_from_repeatmasker(
                repeatmasker_out, record.id, 0, len(seq)
            )
        rows.append(row)
    return rows


FEATURE_COLUMNS = [
    "length",
    "gc_content",
    "female_absent_prop",
    "male_depth",
    "repeat_frac",
]
