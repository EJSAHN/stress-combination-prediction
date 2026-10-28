#!/usr/bin/env python3
"""Freeze fixed four-way expression predictions without parsing target columns.

Only 45 lower-order/discovery library columns are parsed. This module never
opens target numerical columns, fits on them, or scores a prediction against
them. Evaluation is a separate explicitly invoked command.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

CODE2_SHA256 = "4a90353caee5ea2e7a537130fb2fb4425a7776e3debcb8f551e32b9c9d42ef7a"
BASE_MODELS = ["control_mean", "single_effect_sum", "main_effects", "pairwise",
               "ridge_pairwise", "main_clipped", "ridge_pairwise_clipped",
               "rank1_main", "rank1_main_score_clipped", "rank2_main",
               "rank2_main_score_clipped"]
EXTRA_MODELS = ["all_triples_mean", "quiet_triples_mean", "quiet_residual_rank2",
                "quiet_common_rank2", "cubic_factorial"]
MODELS = BASE_MODELS + EXTRA_MODELS
FACTORS = ["salt", "paraquat", "high_light", "heat"]
QUIET = ["PQ-HL-HS", "S-PQ-HL"]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load_reference(path):
    require(sha256(path) == CODE2_SHA256, "CODE2 reference hash mismatch")
    spec = importlib.util.spec_from_file_location("pinned_code2_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def cubic_prediction(means):
    """Degree<=3 factorial extension: triples - pairs + singles - control."""
    return sum(((-1) ** (3 - sum(key))) * value for key, value in means.items())


def synthetic_checks(reference):
    base = reference.synthetic_checks()
    cells = np.array(list(itertools.product((0, 1), repeat=4)))
    columns = [np.ones(len(cells))]
    for degree in (1, 2, 3):
        columns.extend(np.prod(cells[:, indices], axis=1)
                       for indices in itertools.combinations(range(4), degree))
    design = np.column_stack(columns)
    values = design @ (np.arange(15 * 7).reshape(15, 7) / 13)
    means = {tuple(row): value for row, value in zip(cells, values) if sum(row) < 4}
    error = float(np.max(np.abs(cubic_prediction(means) - values[-1])))
    require(error < 1e-10, "Synthetic cubic extension failed")
    basis = np.linalg.qr(np.arange(28).reshape(7, 4) + np.eye(7, 4))[0].T[:2]
    effect = np.arange(7, dtype=float)
    common = (effect @ basis.T) @ basis
    residual = effect - common
    require(np.allclose(common + residual, effect), "Projection identity failed")
    require(np.allclose(np.dot(common, residual), 0), "Projection orthogonality failed")
    return {"code2_checks": base, "cubic_max_error": error,
            "common_residual_dot": float(np.dot(common, residual))}


def validate_protocol(protocol, paths):
    for path in paths:
        require(path.name in protocol["input_hashes"], f"Unpinned input: {path.name}")
        require(sha256(path) == protocol["input_hashes"][path.name],
                f"Input hash mismatch: {path.name}")
    require(protocol["factors"] == FACTORS, "Factor order changed")
    require(protocol["filter_training_mean_count_min"] == 10, "Primary filter changed")
    require(protocol["selected_ridge_alpha"] == 0.01, "CODE2 alpha changed")
    require(protocol["quiet_triples"] == QUIET, "Discovery-defined quiet cells changed")
    require(protocol["models"] == MODELS, "Frozen model roster/order changed")


def model_metadata():
    details = []
    for name in MODELS:
        family = "lower_order_CODE2"
        source = "33 libraries in 11 conditions of order<=2"
        if name in EXTRA_MODELS:
            family = "discovery_informed_extrapolation"
            source = "45 libraries in 15 already-unlocked conditions"
        if name in ("quiet_triples_mean", "quiet_common_rank2", "quiet_residual_rank2"):
            source = "Two quiet triples chosen after CODE2 discovery; lower-order rank2 basis"
        details.append({"model": name, "family": family, "source": source,
                        "outcome_target_fit": False})
    return details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("data", "design", "gene-panel", "code2-module", "protocol", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    reference_path = args.code2_module or Path(__file__).with_name("reference_code2_models.py")
    reference = load_reference(reference_path)
    tests = synthetic_checks(reference)
    if args.self_test:
        print(json.dumps({"status": "PASS", **tests}, indent=2))
        return
    require(all(getattr(args, name) is not None for name in
                ("data", "design", "gene_panel", "protocol", "output")), "Missing arguments")
    protocol = json.loads(args.protocol.read_text())
    inputs = [args.data, args.design, args.gene_panel, reference_path]
    validate_protocol(protocol, inputs)
    protocol_hash = sha256(args.protocol)
    input_hashes = {path.name: sha256(path) for path in inputs}
    code_paths = [Path(__file__), reference_path,
                  Path(__file__).with_name("evaluate_frozen_predictions.py")]
    require(all(path.is_file() for path in code_paths), "Evaluator must exist before code freeze")
    code_hashes = {path.name: sha256(path) for path in code_paths}
    require(not args.output.exists() or not any(args.output.iterdir()),
            "Freeze output must be absent or empty; immutable predictions are never overwritten")
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(args.design)
    require(len(metadata) == 48 and metadata.matrix_column.is_unique, "Invalid sample inventory")
    require(metadata.groupby(FACTORS).size().eq(3).all(), "Expected three pools per factorial cell")
    require(metadata[FACTORS].sum(axis=1).eq(metadata.order).all(), "Invalid factorial order")
    used = metadata.loc[metadata.order < 4].copy().reset_index(drop=True)
    target = metadata.loc[metadata.order == 4].copy().reset_index(drop=True)
    training = used.order.to_numpy() <= 2
    require((len(used), training.sum(), len(target)) == (45, 33, 3), "Invalid frozen split")
    # NEVER PARSE FOUR-WAY NUMERICAL COLUMNS HERE.
    table = pd.read_csv(args.data, sep="\t", usecols=["genes"] + used.matrix_column.tolist())
    require(not set(target.matrix_column).intersection(table.columns), "Target columns were parsed")
    require(table.genes.is_unique, "Duplicate genes")
    counts = table[used.matrix_column].to_numpy(dtype=float).T
    require(np.isfinite(counts).all() and (counts >= 0).all(), "Invalid normalized values")
    keep = counts[training].mean(axis=0) >= 10
    genes = table.loc[keep, "genes"].to_numpy(dtype=str)
    expression = np.log2(1 + counts[:, keep])
    gene_panel = pd.read_csv(args.gene_panel)
    require(np.array_equal(genes, gene_panel.gene.to_numpy(dtype=str)), "CODE2 gene panel differs")
    responsive = gene_panel.training_responsive.to_numpy(dtype=bool)
    require(len(genes) == 21495 and responsive.sum() == 9695, "Unexpected frozen panel sizes")
    factors = used[FACTORS].to_numpy(dtype=int)
    predictions, state = reference.predict_all(expression[training], factors[training],
                                               np.ones((1, 4), dtype=int), 0.01)
    vectors = {name: predictions[name][0] for name in BASE_MODELS}
    means = {key: expression[indices].mean(axis=0)
             for key, indices in reference.cell_indices(factors).items()}
    named_means = {condition: expression[rows.index].mean(axis=0)
                   for condition, rows in used.groupby("condition")}
    control = state["control"]
    basis = state["basis"][:2]
    quiet_full = np.mean([named_means[name] for name in QUIET], axis=0)
    effect = quiet_full - control
    common = (effect @ basis.T) @ basis
    vectors.update({"all_triples_mean": np.mean([value for key, value in means.items()
                                                if sum(key) == 3], axis=0),
                    "quiet_triples_mean": quiet_full,
                    "quiet_common_rank2": control + common,
                    "quiet_residual_rank2": control + effect - common,
                    "cubic_factorial": cubic_prediction(means)})
    require(np.allclose(vectors["quiet_triples_mean"], vectors["quiet_common_rank2"] +
                        vectors["quiet_residual_rank2"] - control), "Quiet identity failed")
    matrix = np.stack([vectors[name] for name in MODELS])
    require(matrix.shape == (16, 21495) and np.isfinite(matrix).all(), "Invalid prediction matrix")
    np.savez_compressed(args.output / "frozen_predictions.npz", genes=genes,
                        model_names=np.array(MODELS), predictions=matrix,
                        control_log_mean=control, training_basis_rank2=basis,
                        training_responsive=responsive,
                        target_columns=target.matrix_column.to_numpy(dtype=str))
    used.to_csv(args.output / "unlocked_sample_design.csv", index=False)
    pd.DataFrame(model_metadata()).to_csv(args.output / "model_metadata.csv", index=False)
    pd.DataFrame({"model": MODELS, "predicted_rms_distance_from_control":
                  np.sqrt(np.mean((matrix - control) ** 2, axis=1)),
                  "predicted_min_log_expression": matrix.min(axis=1),
                  "predicted_max_log_expression": matrix.max(axis=1)}).to_csv(
                      args.output / "prediction_summary.csv", index=False)
    require(sha256(args.protocol) == protocol_hash, "Protocol changed during freezing")
    require(all(sha256(path) == input_hashes[path.name] for path in inputs),
            "Input changed during freezing")
    require(all(sha256(path) == code_hashes[path.name] for path in code_paths),
            "Code changed during freezing")
    manifest = {"status": "PREDICTIONS_FROZEN_NO_TARGET_VALUES_PARSED",
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "protocol_filename": args.protocol.name, "protocol_sha256": protocol_hash,
                "input_hashes": input_hashes,
                "code_hashes": code_hashes,
                "output_hashes": {path.name: sha256(path) for path in args.output.iterdir()
                                  if path.is_file()},
                "unlocked_columns": used.matrix_column.tolist(),
                "target_columns_not_parsed": target.matrix_column.tolist(),
                "models": MODELS, "n_genes": len(genes),
                "n_responsive_genes": int(responsive.sum()),
                "synthetic_checks": tests,
                "caveat": "Internal outcome holdout; source normalization may use all48 libraries."
                " Quiet comparators were selected from explored triples; no independent validation claim."}
    json_write(args.output / "PREDICTIONS_FROZEN.json", manifest)
    print(json.dumps({"status": manifest["status"], "models": 16, "genes": len(genes),
                      "npz_sha256": manifest["output_hashes"]["frozen_predictions.npz"]}))


if __name__ == "__main__":
    main()
