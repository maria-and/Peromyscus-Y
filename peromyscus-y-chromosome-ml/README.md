# Y-DISCOVER-ML: Machine-Learning-Guided Y Chromosome Discovery for Non-Model Genomes

A Snakemake pipeline + trained classifier for finding Y-chromosome-derived scaffolds in a
*de novo* genome assembly of a species with no existing Y reference — built while assembling
the Y chromosome of an F1 hybrid deer mouse (*Peromyscus maniculatus bairdii* x *P. m. nubiterrae*)
for which no Y-bearing reference genome exists in public databases.

## Why this exists

Identifying Y-chromosome sequence in a new genome assembly is a classic non-model-organism
problem: Y is small, repeat-dense, and (unlike autosomes) has no diploid partner to compare
against. The standard approach ([DiscoverY](https://github.com/makovalab-psu/DiscoverY), Rhie
et al. 2019) compares a male assembly's contigs against female sequencing reads — real Y
sequence should be almost entirely **absent** from a female's reads, while autosomal
sequence should be almost entirely **present**. DiscoverY itself uses a hand-set threshold
on this signal.

This project replaces that manual threshold with a **trained classifier**, and adds the
k-mer/depth signal together with repeat content, GC content, and (optionally) homology score
as features — then validates the whole thing end-to-end, because the manual-threshold version
of this analysis is exactly what produced a [596Mb of misclassified autosomal sequence in an
earlier pass of this same project](docs/postmortem_misclassification.md) (see write-up below).
That failure is the actual motivating case study for this repo.

## Pipeline overview

```
raw reads (HiFi + Hi-C/Omni-C)
        │
        ▼
  hifiasm haplotype-resolved assembly  ──▶  per-haplotype scaffolds
        │
        ▼
  feature extraction per contig/window:
    - female-read k-mer absence proportion (Bloom filter)
    - male read depth
    - GC content
    - repeat content (RepeatMasker)
    - [optional] homology score to known Y genes
        │
        ▼
  trained gradient-boosted classifier  ──▶  Y-contig probability per contig
        │
        ▼
  homology + depth cross-validation  ──▶  confirmed Y scaffold set
        │
        ▼
  repeat-masking → BRAKER annotation (RNA-seq + protein evidence)
```

## What makes the classifier different from a manual threshold

Rather than hand-picking a cutoff on one feature (what the earlier failed run effectively did,
and got backwards), `workflow/scripts/train_classifier.py` trains a gradient-boosted tree
classifier on **labeled** contigs pulled from a reference genome with a known, resolved Y
chromosome (e.g. mouse GRCm39: `chrY` vs. autosomes), using the same feature set computed the
same way. The idea is straightforward transfer: the *relationship* between these features and
"is this Y" should generalize across rodent genomes even though the sequence itself doesn't.
The model is evaluated with stratified cross-validation (precision/recall/ROC-AUC, since classes
are wildly imbalanced — Y is a tiny fraction of any genome) before ever touching the actual
target assembly.

## Repo layout

```
config/config.yaml           sample accessions, paths, model hyperparameters
Snakefile                    top-level pipeline entrypoint
workflow/rules/*.smk          one rule file per pipeline stage
workflow/scripts/
  compute_features.py        contig/window feature extraction (reusable, unit-tested)
  build_training_labels.py   pulls labeled training contigs from a reference genome
  train_classifier.py        trains + evaluates the classifier, saves the model + metrics
  classify_contigs.py        applies the trained model to a target assembly's features
  plot_results.py            QC plots: score distributions, ROC curve, confusion matrix
tests/                        pytest unit tests for feature extraction
docs/postmortem_misclassification.md   write-up of the earlier failed manual-threshold run
.github/workflows/ci.yml      lint + test on every push
```

## Running it

```bash
conda env create -f environment.yml
conda activate y-discover-ml

# 1. Build a labeled training set from a reference with a known Y (example: mouse)
python workflow/scripts/build_training_labels.py --reference mouse_GRCm39.fa --out training_labels.csv

# 2. Train + evaluate the classifier
python workflow/scripts/train_classifier.py --labels training_labels.csv --out-model model.pkl

# 3. Run the full pipeline on a target assembly (edit config/config.yaml first)
snakemake --cores 8 --use-conda
```

Or run end-to-end on a compute cluster via Snakemake's Slurm executor:
```bash
snakemake --cores 32 --executor slurm --use-conda
```

## Status

Pipeline structure, feature extraction, and the classifier training/evaluation code are
implemented and unit-tested. Applying this to the actual target assembly (the bairdii x
nubiterrae F1 hybrid) is in progress on a separate HPC cluster — see `docs/postmortem_misclassification.md`
for why the first attempt at this (using a hand-set threshold instead of this classifier) failed.

## License

MIT — see `LICENSE`.
