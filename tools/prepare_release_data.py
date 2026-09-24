"""Export numerical inputs from the original audit, without contributor names.

This maintainer tool prepares a release. Users reproduce the analysis with the
already bundled files and do not need the original workspace.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("pilot", "database", "catpred-variants", "catpred-summary", "out"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        parser.error("Output directory must be new or empty")
    paths = {
        "audit": args.pilot / "01_audit/row_audit.csv",
        "measurements": args.pilot / "02_consensus_dataset/row_deltas.csv",
        "consensus": args.pilot / "02_consensus_dataset/variant_consensus.csv",
        "predictions": args.pilot / "04_baselines/baseline_predictions.csv",
        "reference": args.pilot / "01_audit/wt_reference_audit.json",
        "database": args.database,
        "catpred_variants": args.catpred_variants,
        "catpred_summary": args.catpred_summary,
    }
    audit = read(paths["audit"])
    groups = {}
    for row in audit:
        if row["context_status"] == "PRESENT" and row["context_key"] not in groups:
            groups[row["context_key"]] = f"group_{len(groups) + 1:04d}"

    def group(value):
        if value.startswith("UNVERIFIED::"):
            return value
        return groups.get(value, "UNKNOWN")

    audit_fields = "record_id variant_norm variant_type context_status context_key expressed_yes km_value kcat_value reported_efficiency efficiency_source derived_efficiency reference_km_value reference_kcat_value reported_reference_efficiency reference_efficiency_source derived_reference_efficiency unit_status row_audit_status quarantine_reason_codes flag_reason_codes sequence_validation_status expected_sequence_sha256".split()
    audit_export = [{k: group(row[k]) if k == "context_key" else row[k]
                     for k in audit_fields} for row in audit]
    measurements = []
    for row in read(paths["measurements"]):
        measurements.append({k: group(v) if k in {"context_key", "phase1_context_key"} else v
                             for k, v in row.items() if k != "source_row_number"})
    consensus_fields = "mutation_id position n_rows n_contexts n_verified_contexts quality_tier quality_reason_codes target_delta_log10_eff n_raw_rows source_record_ids target_source_record_ids quarantined_record_ids".split()
    consensus = [{k: row[k] for k in consensus_fields} for row in read(paths["consensus"])]
    variant_results = {r["substitution"]: r for r in read(paths["catpred_variants"])}
    wt_prediction = json.loads(paths["catpred_summary"].read_text())["catpred_wt_median_linear_ratio"]
    predictions = []
    for row in read(paths["predictions"]):
        if row["model"] == "raw_catpred" and row["evaluation"] == "position_holdout":
            result = variant_results[row["mutation_id"]]
            predictions.append({k: row[k] for k in ("mutation_id", "evaluation", "model", "prediction", "observed_delta_log10_eff")}
                               | {"catpred_efficiency_ratio": result["pred_kcat_over_Km"],
                                  "catpred_wt_efficiency_ratio": wt_prediction,
                                  "catpred_input_rows": result["n_rows"]})
    reference_source = json.loads(paths["reference"].read_text())
    reference = {k: reference_source[k] for k in (
        "protein_sequence", "protein_length_aa", "protein_sha256",
        "coding_region_start_1_based", "coding_region_end_1_based_inclusive",
        "coding_region_decision", "trailing_dna", "translation_table")}
    args.out.mkdir(parents=True, exist_ok=True)
    files = {}
    for name, rows in {"audit": audit_export, "measurements": measurements,
                       "consensus": consensus, "predictions": predictions}.items():
        path = args.out / (name + ".csv")
        write(path, rows)
        files[path.name] = {"sha256": digest(path), "rows": len(rows)}
    reference_path = args.out / "reference.json"
    reference_path.write_text(json.dumps(reference, indent=2) + "\n")
    files[reference_path.name] = {"sha256": digest(reference_path)}
    manifest = {
        "description": "Numerical BglB/pNPG analysis inputs derived from public D2D measurements and CatPred results.",
        "data_source": "https://d2d-cure.vercel.app/database/characterization_data/BglB",
        "public_endpoint": "https://d2d-cure.vercel.app/api/getCharacterizationData?enzyme=BglB",
        "snapshot_download_date": None,
        "snapshot_file_modification_date": "2026-02-03",
        "audit_freeze_month": "2026-07",
        "public_source_checked": "2026-09-24",
        "catpred_export_date": "2026-05-15",
        "catpred_code_revision": None,
        "catpred_columns": ["log10kcat_max", "log10km_mean"],
        "catpred_calculation": "10 ** (log10kcat_max - log10km_mean), divided by the WT median; take the median per substitution and then log10.",
        "grouping": "Contributor labels replaced by sequential IDs in source-row order; unknown contributors stay unknown. No name mapping is distributed.",
        "scope": "Reproduces the comparison from audited inputs. Does not rerun CatPred inference or recreate the original raw-data audit.",
        "source_files": {k: {"name": p.name, "sha256": digest(p)} for k, p in paths.items()},
        "files": files,
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({name: info.get("rows") for name, info in files.items()}, indent=2))


if __name__ == "__main__":
    main()
