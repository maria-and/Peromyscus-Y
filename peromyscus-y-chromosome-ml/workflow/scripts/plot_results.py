"""QC plots for the classifier and the final Y-candidate calls.

Usage:
    python plot_results.py --metrics model.pkl.metrics.json --classified y_candidates.csv --out-dir figures/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_feature_importance(metrics: dict, out_path: Path) -> None:
    importance = metrics.get("feature_importance", {})
    if not importance:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    names, values = zip(*sorted(importance.items(), key=lambda kv: kv[1]))
    ax.barh(names, values)
    ax.set_xlabel("Feature importance")
    ax.set_title("Classifier feature importance")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_score_distribution(classified: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(classified["y_probability"], bins=40)
    ax.axvline(0.5, color="red", linestyle="--", label="decision threshold")
    ax.set_xlabel("Predicted Y probability")
    ax.set_ylabel("Number of contigs")
    ax.set_title("Y-probability distribution across target assembly contigs")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_size_vs_score(classified: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = classified["predicted_y"].map({True: "crimson", False: "steelblue"})
    ax.scatter(
        classified["length"], classified["y_probability"], c=colors, alpha=0.6, s=15
    )
    ax.set_xscale("log")
    ax.set_xlabel("Contig length (bp, log scale)")
    ax.set_ylabel("Predicted Y probability")
    ax.set_title(
        "Contig size vs. Y probability\n(the earlier manual-threshold run picked by size, not score)"
    )
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", default=None)
    parser.add_argument("--classified", default=None)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.metrics:
        metrics = json.loads(Path(args.metrics).read_text())
        plot_feature_importance(metrics, out_dir / "feature_importance.png")

    if args.classified:
        classified = pd.read_csv(args.classified)
        plot_score_distribution(classified, out_dir / "score_distribution.png")
        plot_size_vs_score(classified, out_dir / "size_vs_score.png")


if __name__ == "__main__":
    main()
