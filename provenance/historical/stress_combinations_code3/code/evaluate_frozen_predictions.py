#!/usr/bin/env python3
"""Evaluate frozen predictions only after explicit authorization to open target.

This command parses only the three four-way numerical columns. Models, panels,
basis and predictions must already be frozen. It never refits or selects models.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from freeze_predictions import sha256, require, json_write


def metric_rows(predictions, observed, names, columns, genes, responsive, basis):
    rows = []
    panels = {"all_filtered": np.ones(len(genes), dtype=bool),
              "training_responsive": responsive}
    require(np.allclose(basis @ basis.T, np.eye(2), atol=1e-10), "Nonorthogonal basis")
    for model_index, model in enumerate(names):
        error = predictions[model_index] - observed
        common = (error @ basis.T) @ basis
        residual = error - common
        require(np.allclose(np.mean(error**2, axis=1), np.mean(common**2, axis=1) +
                            np.mean(residual**2, axis=1), atol=1e-10),
                "MSE orthogonal decomposition failed")
        for panel, mask in panels.items():
            components = {"expression": error[:, mask]}
            if panel == "all_filtered":
                components.update({"common_rank2": common, "residual_rank2": residual})
            for component, differences in components.items():
                for sample_index, column in enumerate(columns):
                    vector = differences[sample_index]
                    rows.append({"model": str(model), "source_column": str(column),
                                 "panel": panel, "component": component,
                                 "n_genes": int(mask.sum()), "mse": float(np.mean(vector**2)),
                                 "rmse": float(np.sqrt(np.mean(vector**2))),
                                 "mae": float(np.mean(np.abs(vector))),
                                 "role": "internal_single_fourway_condition_3_pooled_libraries"})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("data", "design", "protocol", "freeze-dir", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--authorize-target-read", action="store_true",
                        help="Required explicit opening gate; run only after protocol/predictions frozen")
    args = parser.parse_args()
    require(args.authorize_target_read, "Target remains locked; explicit opening flag is required")
    manifest_path = args.freeze_dir / "PREDICTIONS_FROZEN.json"
    manifest = json.loads(manifest_path.read_text())
    require(manifest["status"] == "PREDICTIONS_FROZEN_NO_TARGET_VALUES_PARSED", "Invalid freeze status")
    require(sha256(args.protocol) == manifest["protocol_sha256"], "Protocol changed after freeze")
    for path in (args.data, args.design):
        require(sha256(path) == manifest["input_hashes"][path.name], f"Changed input: {path.name}")
    for name, expected in manifest["code_hashes"].items():
        require(sha256(Path(__file__).with_name(name)) == expected, f"Changed code: {name}")
    for name, expected in manifest["output_hashes"].items():
        require(sha256(args.freeze_dir / name) == expected, f"Changed frozen output: {name}")
    require(not args.output.exists() or not any(args.output.iterdir()),
            "Evaluation output must be absent or empty; earlier results are never overwritten")
    args.output.mkdir(parents=True, exist_ok=True)
    frozen_path = args.freeze_dir / "frozen_predictions.npz"
    with np.load(frozen_path, allow_pickle=False) as frozen:
        genes = frozen["genes"].copy()
        names = frozen["model_names"].copy()
        predictions = frozen["predictions"].copy()
        control = frozen["control_log_mean"].copy()
        basis = frozen["training_basis_rank2"].copy()
        responsive = frozen["training_responsive"].copy()
        columns = frozen["target_columns"].tolist()
    require(columns == manifest["target_columns_not_parsed"], "Target inventory changed")
    metadata = pd.read_csv(args.design)
    target = metadata.loc[metadata.order == 4].copy().reset_index(drop=True)
    require(target.matrix_column.tolist() == columns and len(columns) == 3, "Invalid target inventory")
    # THIS IS THE SOLE AUTHORIZED TARGET-VALUE READ IN THE MODEL PIPELINE.
    table = pd.read_csv(args.data, sep="\t", usecols=["genes"] + columns)
    require(table.genes.is_unique, "Duplicate target gene identifiers")
    require(set(table.columns) == set(["genes"] + columns), "Unexpected parsed target columns")
    table = table.set_index("genes")
    require(set(genes).issubset(table.index), "Frozen genes absent from target")
    counts = table.loc[genes, columns].to_numpy(dtype=float).T
    require(np.isfinite(counts).all() and (counts >= 0).all(), "Invalid target values")
    observed = np.log2(1 + counts)
    scores = metric_rows(predictions, observed, names, columns, genes, responsive, basis)
    scores.to_csv(args.output / "target_library_metrics.csv", index=False)
    group_keys = ["model", "panel", "component", "n_genes"]
    summary = scores.groupby(group_keys, sort=False).agg(
        n_pooled_libraries=("source_column", "size"), mse_mean=("mse", "mean"),
        mse_min=("mse", "min"), mse_max=("mse", "max"), mae_mean=("mae", "mean"),
        mae_min=("mae", "min"), mae_max=("mae", "max")).reset_index()
    summary.to_csv(args.output / "target_model_comparison.csv", index=False)
    gains = []
    match = ["source_column", "panel", "component"]
    for baseline in ("control_mean", "quiet_common_rank2", "quiet_triples_mean"):
        base = scores.loc[scores.model == baseline, match + ["mse", "mae"]].rename(
            columns={"mse": "baseline_mse", "mae": "baseline_mae"})
        paired = scores.merge(base, on=match, validate="many_to_one")
        paired["baseline"] = baseline
        paired["mse_gain"] = paired.baseline_mse - paired.mse
        paired["mae_gain"] = paired.baseline_mae - paired.mae
        gains.append(paired[["model", "baseline"] + match + ["mse_gain", "mae_gain"]])
    gains = pd.concat(gains, ignore_index=True)
    gains.to_csv(args.output / "paired_library_gains.csv", index=False)
    gains.groupby(["model", "baseline", "panel", "component"], sort=False).agg(
        mse_gain_mean=("mse_gain", "mean"), mse_gain_min=("mse_gain", "min"),
        mse_gain_max=("mse_gain", "max"), mae_gain_mean=("mae_gain", "mean"),
        mae_gain_min=("mae_gain", "min"), mae_gain_max=("mae_gain", "max")
    ).reset_index().to_csv(args.output / "paired_gain_summary.csv", index=False)
    effect = observed - control
    common = (effect @ basis.T) @ basis
    residual = effect - common
    pd.DataFrame({"source_column": columns, "n_genes": len(genes),
                  "total_effect_mse": np.mean(effect**2, axis=1),
                  "common_effect_mse": np.mean(common**2, axis=1),
                  "residual_effect_mse": np.mean(residual**2, axis=1),
                  "total_effect_mae": np.mean(np.abs(effect), axis=1),
                  "common_effect_mae": np.mean(np.abs(common), axis=1),
                  "residual_effect_mae": np.mean(np.abs(residual), axis=1)}).to_csv(
                      args.output / "target_state_geometry.csv", index=False)
    np.savez_compressed(args.output / "evaluated_target.npz", genes=genes,
                        target_columns=np.array(columns), observed_log_expression=observed)
    target.to_csv(args.output / "target_sample_design.csv", index=False)
    result = {"status": "FROZEN_INTERNAL_HOLDOUT_EVALUATED", "n_target_pools": 3,
              "n_target_conditions": 1, "n_models": len(names), "n_genes": len(genes),
              "n_responsive_genes": int(responsive.sum()), "p_values_computed": False,
              "protocol_sha256": sha256(args.protocol), "freeze_manifest_sha256": sha256(manifest_path),
              "frozen_predictions_sha256": sha256(frozen_path),
              "code_hashes": manifest["code_hashes"], "input_hashes": manifest["input_hashes"],
              "outputs": {path.name: sha256(path) for path in args.output.iterdir() if path.is_file()},
              "metric_notes": "MSE/MAE are per pooled library. Positive paired gain means lower error."
              " Common/residual MSE is additive under the fixed orthogonal full-gene basis; MAE is not additive."
              " Responsive genes receive ordinary expression metrics only. No genes-as-replicates inference.",
              "limitations": ["One held-out factorial cell with three pooled libraries, not external validation",
                              "Quiet triple comparators selected after viewing discovery triples",
                              "Author normalization may have used all48 libraries",
                              "Published qualitative four-way findings already known",
                              "No multiplicity-adjusted model-selection inference or P-values"]}
    json_write(args.output / "RUN_SUMMARY.json", result)
    print(json.dumps({key: result[key] for key in ("status", "n_target_pools", "n_models", "n_genes")}))


if __name__ == "__main__":
    main()
