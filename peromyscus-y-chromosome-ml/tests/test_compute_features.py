import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))

from compute_features import (  # noqa: E402
    gc_content,
    kmerize,
    FemaleKmerSet,
    female_absent_proportion,
    windows,
)


def test_gc_content_simple():
    assert gc_content("GGCC") == 1.0
    assert gc_content("AATT") == 0.0
    assert gc_content("ATGC") == 0.5


def test_gc_content_empty_sequence_does_not_divide_by_zero():
    assert gc_content("") == 0.0
    assert gc_content("NNNN") == 0.0


def test_kmerize_basic():
    kmers = kmerize("ACGTACGT", k=4)
    assert "ACGT" in kmers
    assert "CGTA" in kmers
    # 5 positions (8 - 4 + 1), but position 0 and position 4 are both "ACGT" -> 4 unique
    assert len(kmers) == 4


def test_kmerize_excludes_kmers_containing_n():
    kmers = kmerize("ACGNACGT", k=4)
    assert all("N" not in km for km in kmers)


def test_kmerize_short_sequence_returns_empty():
    assert kmerize("AC", k=4) == set()


def test_female_absent_proportion_fully_present_sequence_is_zero():
    female = FemaleKmerSet(k=4)
    female.add_sequence("ACGTACGTACGT")
    # every 4-mer of "ACGT" repeats is already in the female set
    assert female_absent_proportion("ACGTACGT", female) == 0.0


def test_female_absent_proportion_fully_novel_sequence_is_one():
    female = FemaleKmerSet(k=4)
    female.add_sequence("AAAAAAAA")
    assert female_absent_proportion("CCCCGGGG", female) == 1.0


def test_female_absent_proportion_empty_set_treated_as_all_absent():
    female = FemaleKmerSet(k=4)
    assert female_absent_proportion("ACGTACGT", female) == 1.0


def test_windows_covers_full_length_without_gaps_or_overlap():
    spans = list(windows(25, 10))
    assert spans == [(0, 10), (10, 20), (20, 25)]


def test_windows_exact_multiple():
    spans = list(windows(20, 10))
    assert spans == [(0, 10), (10, 20)]
