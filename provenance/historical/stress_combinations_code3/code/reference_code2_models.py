#!/usr/bin/env python3
"""Compare fixed bounded and common-state nulls without reading four-way outcomes.

All responses are continuous log2(1 + author-normalized count). This is an
exploratory predictive diagnostic, not count inference or a mechanistic model.
"""

import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def design_matrix(factors, include_pairs=True):
    factors = np.asarray(factors, dtype=float)
    columns = [np.ones(len(factors))]
    columns.extend(factors[:, index] for index in range(factors.shape[1]))
    if include_pairs:
        columns.extend(
            factors[:, left] * factors[:, right]
            for left, right in itertools.combinations(range(factors.shape[1]), 2)
        )
    return np.column_stack(columns)


def prediction_weights(training_design, target_design, alpha=0.0):
    """Penalize pair coefficients only; retain intercept and four main effects."""
    penalty = np.zeros(training_design.shape[1])
    penalty[5:] = alpha
    return target_design @ np.linalg.solve(
        training_design.T @ training_design + np.diag(penalty),
        training_design.T,
    )


def cell_indices(factors):
    keys = [tuple(int(value) for value in row) for row in factors]
    return {
        key: np.flatnonzero([candidate == key for candidate in keys])
        for key in sorted(set(keys))
    }


def training_state(training_values, training_factors):
    """Control-anchored basis and bounds from lower-order conditions only."""
    groups = cell_indices(training_factors)
    keys = list(groups)
    means = np.stack([training_values[groups[key]].mean(axis=0) for key in keys])
    control = training_values[groups[(0, 0, 0, 0)]].mean(axis=0)
    effects = means - control
    _, singular_values, basis = np.linalg.svd(effects, full_matrices=False)
    # Deterministic orientation for saved scores. Projection is sign invariant.
    for component in range(len(basis)):
        pivot = int(np.argmax(np.abs(basis[component])))
        if basis[component, pivot] < 0:
            basis[component] *= -1
    return {
        "groups": groups,
        "keys": keys,
        "means": means,
        "control": control,
        "minimum": means.min(axis=0),
        "maximum": means.max(axis=0),
        "basis": basis,
        "singular_values": singular_values,
    }


def predict_all(training_values, training_factors, target_factors, alpha):
    state = training_state(training_values, training_factors)
    control = state["control"]
    target_count = len(target_factors)
    main_weights = prediction_weights(
        design_matrix(training_factors, False),
        design_matrix(target_factors, False),
    )
    pair_training = design_matrix(training_factors)
    pair_targets = design_matrix(target_factors)
    predictions = {
        "control_mean": np.repeat(control[None, :], target_count, axis=0),
        "main_effects": main_weights @ training_values,
        "pairwise": prediction_weights(pair_training, pair_targets) @ training_values,
        "ridge_pairwise": prediction_weights(pair_training, pair_targets, alpha)
        @ training_values,
    }
    single_sum = predictions["control_mean"].copy()
    for target_index, target in enumerate(target_factors):
        for factor_index, present in enumerate(target):
            if present:
                single = tuple(int(index == factor_index) for index in range(4))
                single_mean = training_values[state["groups"][single]].mean(axis=0)
                single_sum[target_index] += single_mean - control
    predictions["single_effect_sum"] = single_sum
    for original, bounded in [
        ("main_effects", "main_clipped"),
        ("ridge_pairwise", "ridge_pairwise_clipped"),
    ]:
        predictions[bounded] = np.clip(
            predictions[original], state["minimum"], state["maximum"]
        )
    for rank in (1, 2):
        basis = state["basis"][:rank]
        library_scores = (training_values - control) @ basis.T
        predicted_scores = main_weights @ library_scores
        mean_scores = (state["means"] - control) @ basis.T
        predictions[f"rank{rank}_main"] = control + predicted_scores @ basis
        bounded_scores = np.clip(
            predicted_scores, mean_scores.min(axis=0), mean_scores.max(axis=0)
        )
        predictions[f"rank{rank}_main_score_clipped"] = control + bounded_scores @ basis
    return predictions, state


def synthetic_checks():
    """Check exact factorial recovery, rank reconstruction and prediction bounds."""
    cells = np.array(list(itertools.product((0, 1), repeat=4)))
    training = cells.sum(axis=1) <= 2
    target = cells.sum(axis=1) == 3
    matrix = design_matrix(cells)
    coefficients = np.arange(11 * 8, dtype=float).reshape(11, 8) / 17
    values = matrix @ coefficients
    predictions, state = predict_all(values[training], cells[training], cells[target], 0)
    error = float(np.max(np.abs(predictions["pairwise"] - values[target])))
    require(error < 1e-10, "Synthetic pairwise extrapolation failed")
    require(np.linalg.matrix_rank(matrix[training]) == 11, "Training design rank failed")
    for model in ("main_clipped", "ridge_pairwise_clipped"):
        require(np.all(predictions[model] >= state["minimum"] - 1e-10), "Lower bound failed")
        require(np.all(predictions[model] <= state["maximum"] + 1e-10), "Upper bound failed")
    loading = np.linspace(-1, 1, 8)
    scores = cells @ np.array([0.3, -0.2, 0.4, 0.1])
    rank_one_values = 4.0 + scores[:, None] * loading[None, :]
    rank_predictions, _ = predict_all(
        rank_one_values[training], cells[training], cells[target], 0
    )
    rank_error = float(np.max(np.abs(rank_predictions["rank1_main"] - rank_one_values[target])))
    require(rank_error < 1e-10, "Synthetic rank-one main-effect prediction failed")
    return {"pairwise_max_error": error, "rank_one_max_error": rank_error, "rank": 11}


def score_predictions(predictions, observed, sample_metadata, panels, control):
    rows = []
    for model, prediction in predictions.items():
        require(np.isfinite(prediction).all(), f"Nonfinite predictions: {model}")
        for panel, mask in panels.items():
            differences = prediction[:, mask] - observed[:, mask]
            baseline_mse = ((observed[:, mask] - control[mask]) ** 2).mean(axis=1)
            for index, sample in sample_metadata.iterrows():
                mse = float(np.mean(differences[index] ** 2))
                rows.append(
                    {
                        "condition": sample["condition"],
                        "GSM": sample["GSM"],
                        "source_column": sample["matrix_column"],
                        "model": model,
                        "panel": panel,
                        "n_genes": int(mask.sum()),
                        "mse": mse,
                        "rmse": float(np.sqrt(mse)),
                        "mae": float(np.mean(np.abs(differences[index]))),
                        "control_baseline_mse": float(baseline_mse[index]),
                        "mse_fraction_of_control_baseline": mse / baseline_mse[index],
                        "role": "already_explored_triple_discovery",
                    }
                )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("data", "design", "protocol", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    tests = synthetic_checks()
    if arguments.self_test:
        print(json.dumps({"status": "PASS", **tests}, indent=2))
        return
    require(all(getattr(arguments, name) is not None for name in
                ("data", "design", "protocol", "output")), "Missing path arguments")
    started = time.time()
    protocol = json.loads(arguments.protocol.read_text())
    for path in (arguments.data, arguments.design):
        require(file_sha256(path) == protocol["input_hashes"][path.name],
                f"Input checksum mismatch: {path.name}")
    arguments.output.mkdir(parents=True, exist_ok=True)
    factors = protocol["factors"]
    metadata = pd.read_csv(arguments.design)
    require(metadata.groupby(factors).size().eq(3).all(), "Expected three libraries per cell")
    require(len(metadata[factors].drop_duplicates()) == 16, "Incomplete factorial design")
    require(metadata.matrix_column.is_unique, "Duplicated sample labels")
    require(metadata[factors].sum(axis=1).eq(metadata.order).all(), "Incorrect factor order")
    used = metadata.loc[metadata.order < 4].copy().reset_index(drop=True)
    locked = metadata.loc[metadata.order == 4, "matrix_column"].tolist()
    training_mask = used.order.to_numpy() <= 2
    target_mask = used.order.to_numpy() == 3
    require((training_mask.sum(), target_mask.sum(), len(locked)) == (33, 12, 3), "Incorrect split")
    # Four-way numerical columns are excluded before parsing arrays.
    table = pd.read_csv(arguments.data, sep="\t", usecols=["genes"] + used.matrix_column.tolist())
    require(not set(locked).intersection(table.columns), "Locked columns were parsed")
    require(table.genes.is_unique, "Duplicate gene identifiers")
    counts = table[used.matrix_column].to_numpy(dtype=float).T
    require(np.isfinite(counts).all() and np.all(counts >= 0), "Invalid normalized counts")
    keep = counts[training_mask].mean(axis=0) >= protocol["filter_training_mean_count_min"]
    genes = table.loc[keep, "genes"].to_numpy(dtype=str)
    expression = np.log2(1.0 + counts[:, keep])
    factor_values = used[factors].to_numpy(dtype=int)
    groups = cell_indices(factor_values)
    control_key = (0, 0, 0, 0)
    means = {key: expression[index].mean(axis=0) for key, index in groups.items()}
    variances = {key: expression[index].var(axis=0, ddof=1) for key, index in groups.items()}
    control = means[control_key]
    responsive = np.zeros(len(genes), dtype=bool)
    for key in groups:
        if 0 < sum(key) <= 2:
            effect = means[key] - control
            responsive |= (np.abs(effect) >= 1) & (effect**2 >
                           4 * (variances[key] + variances[control_key]) / 3)
    panels = {"all_filtered": np.ones(len(genes), dtype=bool), "training_responsive": responsive}
    pd.DataFrame({"gene": genes, "training_responsive": responsive,
                  "training_mean_normalized_count": counts[training_mask][:, keep].mean(axis=0)}
                 ).to_csv(arguments.output / "gene_panel.csv", index=False)
    folds = np.full(len(used), -1, dtype=int)
    for key, indices in groups.items():
        if sum(key) <= 2:
            folds[indices] = np.arange(len(indices))
    ridge_rows = []
    pair_design = design_matrix(factor_values)
    for alpha in protocol["ridge_alpha_grid"]:
        for fold in range(3):
            fit = training_mask & (folds != fold)
            validation = training_mask & (folds == fold)
            predicted = prediction_weights(pair_design[fit], pair_design[validation], alpha) @ expression[fit]
            for panel, mask in panels.items():
                ridge_rows.append({"alpha": alpha, "fold": fold, "panel": panel,
                    "mean_library_mse": float(np.mean((predicted[:, mask] - expression[validation][:, mask]) ** 2))})
    ridge_cv = pd.DataFrame(ridge_rows)
    ridge_cv.to_csv(arguments.output / "ridge_training_cv.csv", index=False)
    ranking = (ridge_cv.loc[ridge_cv.panel == "all_filtered"]
               .groupby("alpha", as_index=False).mean_library_mse.mean()
               .sort_values(["mean_library_mse", "alpha"]))
    chosen_alpha = float(ranking.iloc[0].alpha)
    # Known-condition replicate prediction is a denoising diagnostic only.
    model_cv_rows = []
    for fold in range(3):
        fit = training_mask & (folds != fold)
        validation = training_mask & (folds == fold)
        fold_predictions, _ = predict_all(
            expression[fit], factor_values[fit], factor_values[validation], chosen_alpha)
        for model, prediction in fold_predictions.items():
            for panel, mask in panels.items():
                model_cv_rows.append({"fold": fold, "model": model, "panel": panel,
                    "mean_library_mse": float(np.mean((prediction[:, mask] - expression[validation][:, mask]) ** 2)),
                    "role": "training_replicate_denoising_diagnostic_not_outer_validation"})
    pd.DataFrame(model_cv_rows).to_csv(arguments.output / "model_training_cv.csv", index=False)
    predictions, state = predict_all(expression[training_mask], factor_values[training_mask],
                                    factor_values[target_mask], chosen_alpha)
    require(set(predictions) == set(protocol["models"]), "Model set differs from frozen protocol")
    target_metadata = used.loc[target_mask].reset_index(drop=True)
    scores = score_predictions(predictions, expression[target_mask], target_metadata, panels, control)
    scores.to_csv(arguments.output / "discovery_library_scores.csv", index=False)
    comparison = scores.groupby(["condition", "model", "panel"], as_index=False).agg(
        n_libraries=("GSM", "size"), n_genes=("n_genes", "first"),
        mean_mse=("mse", "mean"), mean_mae=("mae", "mean"),
        minimum_library_rmse=("rmse", "min"), maximum_library_rmse=("rmse", "max"),
        baseline_mse=("control_baseline_mse", "mean"))
    comparison["pooled_rmse"] = np.sqrt(comparison.mean_mse)
    comparison["mse_fraction_of_control_baseline"] = comparison.mean_mse / comparison.baseline_mse
    comparison["role"] = "already_explored_triple_discovery"
    comparison.to_csv(arguments.output / "discovery_model_comparison.csv", index=False)
    overall = comparison.groupby(["model", "panel"], as_index=False).agg(
        n_conditions=("condition", "size"), equal_condition_mse=("mean_mse", "mean"),
        equal_condition_mae=("mean_mae", "mean"))
    overall["pooled_rmse"] = np.sqrt(overall.equal_condition_mse)
    overall.to_csv(arguments.output / "model_overview.csv", index=False)
    energy = state["singular_values"] ** 2
    basis_rows = [{"component": index + 1,
                   "training_effect_energy_fraction": float(value / energy.sum()),
                   "cumulative_energy_fraction": float(energy[:index + 1].sum() / energy.sum()),
                   "scale": "unscaled_log2_effects_control_anchored"}
                  for index, value in enumerate(energy)]
    pd.DataFrame(basis_rows).to_csv(arguments.output / "training_basis.csv", index=False)
    # Absolute-state checks distinguish conditional reversal from a new state.
    state_rows = []
    for key, indices in groups.items():
        effect = means[key] - control
        total_energy = float(effect @ effect)
        row = {"condition": used.loc[indices[0], "condition"], "order": sum(key),
               "role": "training" if sum(key) <= 2 else "discovery",
               "rms_distance_from_control": float(np.sqrt(np.mean(effect**2)))}
        for rank in (1, 2):
            basis = state["basis"][:rank]
            residual = effect - (effect @ basis.T) @ basis
            row[f"rank{rank}_absolute_state_residual_rms"] = float(np.sqrt(np.mean(residual**2)))
            row[f"rank{rank}_absolute_state_residual_energy_fraction"] = (
                float(residual @ residual / total_energy) if total_energy > 0 else 0.0)
        state_rows.append(row)
    pd.DataFrame(state_rows).to_csv(arguments.output / "absolute_state_geometry.csv", index=False)
    conditions = target_metadata.condition.drop_duplicates().tolist()
    first_rows = [int(target_metadata.index[target_metadata.condition == name][0]) for name in conditions]
    prediction_archive = {"genes": genes, "conditions": np.array(conditions, dtype=str),
                          "control_log_mean": control, "training_responsive": responsive,
                          "training_basis_rank2": state["basis"][:2]}
    prediction_archive.update({model: value[first_rows] for model, value in predictions.items()})
    np.savez_compressed(arguments.output / "condition_predictions.npz", **prediction_archive)
    used.assign(training_cv_fold=folds).to_csv(arguments.output / "used_sample_design.csv", index=False)
    summary = {
        "status": "EXPLORATORY_BASELINE_COMPARISON_COMPLETE_DATASET_RETAINED",
        "source_genes": len(table), "filtered_genes": len(genes),
        "training_responsive_genes": int(responsive.sum()), "models": protocol["models"],
        "training_libraries": 33, "triple_discovery_libraries": 12,
        "locked_columns": locked, "locked_expression_columns_read": False,
        "selected_ridge_alpha": chosen_alpha, "synthetic_checks": tests,
        "training_rank1_energy_fraction": float(energy[0] / energy.sum()),
        "training_rank2_energy_fraction": float(energy[:2].sum() / energy.sum()),
        "protocol_sha256": file_sha256(arguments.protocol), "code_sha256": file_sha256(__file__),
        "input_sha256": {path.name: file_sha256(path) for path in (arguments.data, arguments.design)},
        "limits": protocol["interpretation"], "runtime_seconds": round(time.time() - started, 3),
    }
    (arguments.output / "RUN_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
