"""Reproduce the real BglB comparison from the numerical inputs in this repository."""

import argparse
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

from .analysis import number, read_csv, run, unique, write_csv


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(left, right, label):
    require(math.isclose(number(left), number(right), rel_tol=1e-12, abs_tol=1e-12),
            f"Does not reproduce: {label}")


def positive(value):
    return bool(value) and number(value) > 0


def efficiency(row, reference=False):
    prefix = "reference_" if reference else ""
    kcat, km = row[prefix + "kcat_value"], row[prefix + "km_value"]
    reported = row["reported_reference_efficiency" if reference else "reported_efficiency"]
    if positive(kcat) and positive(km):
        return number(kcat) / number(km), "kcat_div_KM"
    if positive(reported):
        return number(reported), "reported_kcatKM"
    return None, "missing"


def ids(value):
    return set(filter(None, value.split("|")))


def validate(data):
    """Check fixed inputs, then rebuild efficiencies, WT comparisons and targets.

    Original exclusion and sequence-audit decisions are inputs, not newly
    inferred experimental facts. The checks below verify their consistency.
    """
    data = Path(data)
    manifest = json.loads((data / "manifest.json").read_text())
    expected_files = {"audit.csv", "measurements.csv", "consensus.csv", "predictions.csv", "reference.json"}
    require(set(manifest["files"]) == expected_files, "Unexpected input file list")
    for name, info in manifest["files"].items():
        require(hashlib.sha256((data / name).read_bytes()).hexdigest() == info["sha256"],
                f"Input hash changed: {name}")
    tables = {name: read_csv(data / (name + ".csv"), [])
              for name in ("audit", "measurements", "consensus", "predictions")}
    for name, rows in tables.items():
        require(len(rows) == manifest["files"][name + ".csv"]["rows"], f"Row count changed: {name}")
    audit = unique(tables["audit"], "record_id")
    measurements = unique(tables["measurements"], "record_id")
    consensus = unique(tables["consensus"], "mutation_id")
    predictions = unique(tables["predictions"], "mutation_id")
    reference = json.loads((data / "reference.json").read_text())
    wt = reference["protein_sequence"]
    require(len(wt) == reference["protein_length_aa"] and
            hashlib.sha256(wt.encode()).hexdigest() == reference["protein_sha256"], "WT sequence mismatch")
    observed_eff, reference_eff = {}, {}
    by_mutation = defaultdict(set)
    wt_by_group = defaultdict(set)
    wt_ids = set()
    for identifier, row in audit.items():
        for ref, output in ((False, observed_eff), (True, reference_eff)):
            value, source = efficiency(row, reference=ref)
            source_field = "reference_efficiency_source" if ref else "efficiency_source"
            value_field = "derived_reference_efficiency" if ref else "derived_efficiency"
            require(row[source_field] == source, f"Efficiency source mismatch: {identifier}")
            if value is None:
                require(not row[value_field], f"Unexpected efficiency: {identifier}")
            else:
                close(value, row[value_field], f"efficiency in record {identifier}")
            output[identifier] = value
        mutation = row["variant_norm"]
        sequence = wt
        if row["variant_type"] == "SINGLE":
            position = int(mutation[1:-1]) - 1
            require(0 <= position < len(wt) and wt[position] == mutation[0]
                    and mutation[-1] != wt[position], f"Mutation does not match WT: {mutation}")
            sequence = wt[:position] + mutation[-1] + wt[position + 1:]
            by_mutation[mutation].add(identifier)
        else:
            require(row["variant_type"] == "WT", f"Unexpected variant type: {identifier}")
        require(hashlib.sha256(sequence.encode()).hexdigest() == row["expected_sequence_sha256"],
                f"Expected sequence hash mismatch: {identifier}")
        status = row["row_audit_status"]
        require(status in {"ELIGIBLE", "REFERENCE_ONLY", "QUARANTINE"}, f"Unknown audit status: {identifier}")
        if status == "QUARANTINE":
            require(bool(row["quarantine_reason_codes"]), f"Missing exclusion reason: {identifier}")
        else:
            require(row["expressed_yes"].lower() == "true" and observed_eff[identifier] is not None
                    and row["unit_status"] == "PASS_HEADER_DECLARED"
                    and not row["quarantine_reason_codes"], f"Unusable record was retained: {identifier}")
        if status == "REFERENCE_ONLY":
            require(row["variant_type"] == "WT", f"WT baseline uses mutant: {identifier}")
            wt_ids.add(identifier)
            if row["context_status"] in {"PRESENT", "DOCUMENTED"}:
                wt_by_group[row["context_key"]].add(identifier)
    eligible_ids = {key for key, row in audit.items() if row["row_audit_status"] == "ELIGIBLE"}
    require(eligible_ids == set(measurements), "Eligible audit records do not match measurements")
    for identifier, row in measurements.items():
        source = audit[identifier]
        require(row["mutation_id"] == source["variant_norm"] and source["variant_type"] == "SINGLE",
                f"Measurement mutation mismatch: {identifier}")
        require(row["context_status"] == source["context_status"], f"Contributor status mismatch: {identifier}")
        expected_group = (f"UNVERIFIED::{identifier}" if source["context_status"] == "PLACEHOLDER_UNKNOWN"
                          else source["context_key"])
        require(row["context_key"] == expected_group, f"Contributor grouping mismatch: {identifier}")
        if reference_eff[identifier] is not None:
            baseline_type, baseline_ids = "reference_wt", {identifier}
            baseline = reference_eff[identifier]
        else:
            group_ids = wt_by_group[source["context_key"]]
            baseline_type, baseline_ids = ("context_wt", group_ids) if group_ids else ("global_wt", wt_ids)
            baseline = statistics.median(observed_eff[key] for key in baseline_ids)
        require(row["baseline_type"] == baseline_type and ids(row["baseline_source_record_ids"]) == baseline_ids,
                f"WT source mismatch: {identifier}")
        close(baseline, row["baseline_value"], f"WT baseline for record {identifier}")
        close(observed_eff[identifier], row["efficiency_value"], f"measurement {identifier}")
        close(math.log10(observed_eff[identifier] / baseline), row["delta_log10_eff"], f"effect for record {identifier}")
    require(set(consensus) == set(by_mutation), "Consensus is missing a source mutation")
    for mutation, row in consensus.items():
        source_ids = by_mutation[mutation]
        usable = source_ids & eligible_ids
        require(ids(row["source_record_ids"]) == source_ids and int(row["n_raw_rows"]) == len(source_ids),
                f"Total records do not reconcile: {mutation}")
        require(ids(row["target_source_record_ids"]) == usable and
                ids(row["quarantined_record_ids"]) == source_ids - usable,
                f"Excluded records do not reconcile: {mutation}")
    eligible_mutations = {row["mutation_id"] for row in measurements.values()}
    require(set(predictions) == eligible_mutations, "CatPred coverage does not match eligible mutations")
    for mutation, row in predictions.items():
        require(row["model"] == "raw_catpred" and row["evaluation"] == "position_holdout", "Unexpected prediction model")
        close(row["observed_delta_log10_eff"], consensus[mutation]["target_delta_log10_eff"], f"prediction target {mutation}")
        close(math.log10(number(row["catpred_efficiency_ratio"]) / number(row["catpred_wt_efficiency_ratio"])),
              row["prediction"], f"CatPred normalization {mutation}")
    return {"source_records": len(audit), "usable_wt_records": len(wt_ids),
            "usable_mutant_records": len(measurements), "excluded_source_records": len(audit) - len(wt_ids) - len(measurements),
            "mutations": len(consensus), "catpred_results": len(predictions)}, audit


def reproduce(data, policy, out, plot=False):
    data, out = Path(data), Path(out)
    checks, audit = validate(data)
    summary = run(argparse.Namespace(
        consensus=data / "consensus.csv", measurements=data / "measurements.csv",
        predictions=data / "predictions.csv", policy=policy, out=out, model="raw_catpred",
        dataset_label="D2D BglB measurements and CatPred results", plot=plot))
    records = read_csv(out / "record_selection.csv", [])
    for row in records:
        original = audit[row["record_id"]]
        row["original_audit_status"] = original["row_audit_status"]
        row["original_exclusion_reasons"] = original["quarantine_reason_codes"]
        if row["row_exclusion_reason"] == "excluded_in_earlier_raw_data_audit":
            row["row_exclusion_reason"] = "earlier_audit:" + original["quarantine_reason_codes"]
    write_csv(out / "record_selection.csv", records)
    summary["data_validation"] = checks
    summary["input_manifest_sha256"] = hashlib.sha256((data / "manifest.json").read_bytes()).hexdigest()
    summary["software"]["reproduction_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/bglb"), help="Bundled numerical data directory")
    parser.add_argument("--policy", type=Path, default=Path("selection_policy.json"))
    parser.add_argument("--out", type=Path, default=Path("results"))
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    try:
        summary = reproduce(args.data, args.policy, args.out, args.plot)
    except (ValueError, KeyError, OSError, ImportError) as error:
        parser.exit(2, f"Cannot reproduce comparison: {error}\n")
    print(json.dumps({"selection": summary["selection"], "data_validation": summary["data_validation"]}, indent=2))


if __name__ == "__main__":
    main()
