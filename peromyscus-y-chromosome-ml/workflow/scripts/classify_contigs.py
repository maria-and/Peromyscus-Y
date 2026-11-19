"""
Apply a trained classifier to a target assembly's contigs/windows, producing a
ranked list of Y-candidate contigs with calibrated probabilities instead of a
single hand-set threshold.

Usage:
    python classify_contigs.py \
        --model model.pkl \
        --features target_features.csv \
        --out y_candidates.csv \
        --known-y-genes-gff known_y_hits.gff   # optional cross-validation
"""

from __future__ import annotations

import argparse
import pickle
import sys

import pandas as pd


def load_model(model_path: str):
    with open(model_path, "rb") as fh:
        bundle = pickle.load(fh)
    return bundle["model"], bundle["feature_columns"]


def classify(
    features_df: pd.DataFrame, model, feature_columns: list[str]
) -> pd.DataFrame:
    missing = [c for c in feature_columns if c not in features_df.columns]
    if missing:
        raise ValueError(f"features file is missing required columns: {missing}")

    X = features_df[feature_columns]
    proba = model.predict_proba(X)[:, 1]
    out = features_df.copy()
    out["y_probability"] = proba
    out["predicted_y"] = proba >= 0.5
    return out.sort_values("y_probability", ascending=False)


def cross_validate_with_homology(
    classified: pd.DataFrame, homology_gff: str
) -> pd.DataFrame:
    """Flag which classifier-predicted-Y contigs also have a hit in a homology
    GFF (e.g. miniprot output against known Y genes). Agreement between the
    classifier and independent homology evidence is the real confidence signal
    — see docs/postmortem_misclassification.md for why relying on one signal
    alone isn't enough."""
    hit_contigs = set()
    with open(homology_gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            hit_contigs.add(line.split("\t")[0])
    classified = classified.copy()
    classified["homology_confirmed"] = classified["contig"].isin(hit_contigs)
    return classified


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True)
    parser.add_argument(
        "--features", required=True, help="CSV with a 'contig' column + features"
    )
    parser.add_argument("--out", required=True)
    parser.add_argument("--known-y-genes-gff", default=None)
    args = parser.parse_args()

    model, feature_columns = load_model(args.model)
    features_df = pd.read_csv(args.features)
    classified = classify(features_df, model, feature_columns)

    if args.known_y_genes_gff:
        classified = cross_validate_with_homology(classified, args.known_y_genes_gff)
        disagreements = classified[
            classified["predicted_y"] & ~classified["homology_confirmed"]
        ]
        if len(disagreements) > 0:
            print(
                f"NOTE: {len(disagreements)} contig(s) predicted Y by the classifier have "
                "no homology confirmation yet — worth manual review before trusting them.",
                file=sys.stderr,
            )

    classified.to_csv(args.out, index=False)
    n_predicted = int(classified["predicted_y"].sum())
    print(
        f"Wrote {len(classified)} scored contigs ({n_predicted} predicted Y) to {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
