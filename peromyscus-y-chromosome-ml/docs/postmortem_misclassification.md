# Postmortem: how a manual-threshold DiscoverY run misclassified 609Mb of autosomal sequence as "Y"

This is a real failure from an earlier pass of this project, kept here as the motivating case
study for replacing manual thresholds with a trained, evaluated classifier.

## What happened

A prior run selected 7 large contigs (609.4Mb total, largest single contig 204.3Mb) from a
*de novo* hifiasm assembly as "the Y chromosome," based on a DiscoverY-style per-contig score:
the proportion of each contig's k-mers *absent* from a female short-read dataset, plus male
read depth.

## Why it was wrong

1. **Scale**: 609Mb is chromosome-scale, not Y-scale. A real Y chromosome in this size class of
   rodent should be tens of Mb, not hundreds.
2. **The score itself was backwards for these contigs**: real Y-specific sequence should show a
   *high* proportion of k-mers absent from female data (females have no Y at all). The selected
   contigs showed only 0.15–0.18 — meaning 82–85% of their k-mers *were* found in female data,
   i.e. this sequence is mostly **shared with the female genome**, the opposite of a Y signature.
3. **Genome-wide context made this obvious in hindsight**: across all 440 contigs in that
   assembly, the single highest female-absence score anywhere was 0.21 — nothing approaching the
   0.8–0.9+ expected for real Y sequence. The 7 selected contigs weren't even the highest-scoring
   ones available; the top-scoring contig in the whole assembly wasn't included at all. The
   final selection didn't consistently follow the computed scores — it looks like contig *size*
   drove the pick more than the score did.
4. **Gene content confirmed it independently**: of 190 genes annotated on these "Y" contigs
   (via RNA-seq from male liver, StringTie), none were canonical Y genes (no *Sry*, *Zfy*,
   *Ddx3y*, *Uty*, *Eif2s3y*, *Usp9y*, *Kdm5d*) — instead generic genes like *teneurin-3* and
   *Wnt-3*, consistent with ordinary autosomal content.
5. **A likely secondary bug**: several contigs in the genome-wide score table showed a "depth"
   value of ~13,000,000–13,600,000 — not a real sequencing depth (real values elsewhere were
   40–70x). This looks like a placeholder/error value leaking out of whatever script computed
   depth, suggesting the underlying computation had a real bug independent of the threshold
   question.

## The lesson this repo is built around

A single hand-picked threshold on one feature, applied without:
- checking the score distribution across the *entire* genome (not just the pre-selected
  candidates),
- cross-validating against known gene homology, and
- sanity-checking the absolute scale of the result against prior biological expectation,

...can silently produce a confident-looking but completely wrong answer. This repo's classifier
is trained and evaluated (precision/recall/ROC-AUC under cross-validation, on a labeled
reference genome) specifically so there's a quantified, checkable basis for any Y/not-Y call,
instead of a single eyeballed cutoff.
