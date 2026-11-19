"""
Train and evaluate a gradient-boosted classifier for Y-vs-not-Y window/contig
classification, replacing a hand-set threshold with a model whose precision/
recall/ROC-AUC are actually measured under cross-validation.

Y is a tiny fraction of any genome, so this uses stratified k-fold CV and
reports precision/recall/F1/ROC-AUC (not accuracy, which would be ~99%+ even
for a model that calls everything "not Y").

Usage:
    python train_classifier.py --labels training_labels.csv --out-model model.pkl
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys

import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
)

FEATURE_COLUMNS = ["length", "gc_content", "female_absent_prop", "male_depth"]


def load_training_data(labels_csv: str) -> tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(labels_csv)
    df = df.dropna(subset=FEATURE_COLUMNS + ["label"])
    X = df[FEATURE_COLUMNS]
    y = df["label"].astype(int)
    return X, y


def evaluate_cv(
    X: pd.DataFrame, y: pd.Series, model, cv_folds: int, random_state: int
) -> dict:
    skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    y_pred = cross_val_predict(model, X, y, cv=skf)
    y_proba = cross_val_predict(model, X, y, cv=skf, method="predict_proba")[:, 1]

    tn, fp, fn, tp = confusion_matrix(y, y_pred).ravel()
    return {
        "precision": precision_score(y, y_pred, zero_division=0),
        "recall": recall_score(y, y_pred, zero_division=0),
        "f1": f1_score(y, y_pred, zero_division=0),
        "roc_auc": roc_auc_score(y, y_proba) if len(set(y)) > 1 else float("nan"),
        "confusion_matrix": {
            "tn": int(tn),
            "fp": int(fp),
            "fn": int(fn),
            "tp": int(tp),
        },
        "n_examples": len(y),
        "n_positive": int(y.sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--out-model", required=True)
    parser.add_argument(
        "--out-metrics", default=None, help="Defaults to <out-model>.metrics.json"
    )
    parser.add_argument("--n-estimators", type=int, default=300)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    X, y = load_training_data(args.labels)
    if y.sum() < args.cv_folds or (len(y) - y.sum()) < args.cv_folds:
        print(
            "WARNING: too few positive/negative examples for the requested CV fold count — "
            "reduce --cv-folds or add more labeled data.",
            file=sys.stderr,
        )

    model = GradientBoostingClassifier(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        random_state=args.random_state,
    )

    metrics = evaluate_cv(X, y, model, args.cv_folds, args.random_state)
    print(json.dumps(metrics, indent=2), file=sys.stderr)

    # Fit on all data for the final saved model (CV above is purely for evaluation)
    model.fit(X, y)
    feature_importance = dict(zip(FEATURE_COLUMNS, model.feature_importances_.tolist()))
    metrics["feature_importance"] = feature_importance

    with open(args.out_model, "wb") as fh:
        pickle.dump({"model": model, "feature_columns": FEATURE_COLUMNS}, fh)

    metrics_path = args.out_metrics or f"{args.out_model}.metrics.json"
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)

    print(
        f"Saved model to {args.out_model}, metrics to {metrics_path}", file=sys.stderr
    )


if __name__ == "__main__":
    main()
