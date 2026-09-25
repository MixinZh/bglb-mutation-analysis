"""A small, read-only consumer of previously audited BglB tables.

Selection is calculated before predictions are joined. Original measurements,
quality decisions, and model predictions are never changed by this program.
"""

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def read_csv(path, required):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate column in {Path(path).name}")
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f"Missing columns in {Path(path).name}: {sorted(missing)}")
        rows = list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError(f"Malformed row in {Path(path).name}")
        return rows


def number(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Numeric inputs must be finite")
    return result


def unique(rows, key):
    result = {}
    for row in rows:
        value = row[key]
        if not value or value in result:
            raise ValueError(f"Missing or duplicate {key}: {value!r}")
        result[value] = row
    return result


def load_policy(path):
    policy = json.loads(Path(path).read_text())
    for field in ("minimum_eligible_records", "minimum_contributor_groups"):
        if type(policy[field]) is not int or policy[field] < 4:
            raise ValueError(f"{field} must be an integer of at least four")
    for field in ("effect_threshold_log10", "high_spread_threshold_log10"):
        if number(policy[field]) <= 0:
            raise ValueError(f"{field} must be positive")
    if type(policy["highlight_minimum_records"]) is not int or policy["highlight_minimum_records"] < policy["minimum_eligible_records"]:
        raise ValueError("Highlight count must be an integer at least as large as the record minimum")
    return policy


def record_exclusion(row):
    reasons = []
    if row["baseline_type"] not in {"reference_wt", "context_wt"}:
        reasons.append("no_matched_wt_reference")
    if row["context_status"] not in {"PRESENT", "DOCUMENTED"}:
        reasons.append("unknown_contributor")
    if row["phase1_downweight_reason_codes"]:
        reasons.append("earlier_row_quality_flag:" + row["phase1_downweight_reason_codes"])
    return "|".join(reasons)


def summarize_records(records):
    contexts = defaultdict(list)
    for row in records:
        contexts[row["context_key"]].append(number(row["delta_log10_eff"]))
    medians = [statistics.median(values) for values in contexts.values()]
    return (statistics.median(medians) if medians else None), medians


def select(consensus, measurements, policy):
    """Return the complete selection ledger without accepting model inputs."""
    variants = unique(consensus, "mutation_id")
    unique(measurements, "record_id")
    grouped = defaultdict(list)
    for row in measurements:
        if row["mutation_id"] not in variants:
            raise ValueError("Measurement mutation missing from consensus")
        if not row["context_key"]:
            raise ValueError("Measurement needs a contributor-group key")
        if row["context_status"] not in {"PRESENT", "DOCUMENTED", "PLACEHOLDER_UNKNOWN"}:
            raise ValueError("Unknown context status")
        number(row["delta_log10_eff"])
        grouped[row["mutation_id"]].append(row)
    ledger = []
    for mutation, row in variants.items():
        if not re.fullmatch(r"[A-Z][1-9][0-9]*[A-Z]", mutation):
            raise ValueError(f"Invalid single-substitution label: {mutation}")
        if int(mutation[1:-1]) != int(row["position"]) or mutation[0] == mutation[-1]:
            raise ValueError(f"Mutation label and position disagree: {mutation}")
        if row["quality_tier"] not in {"CORE", "EXTENDED", "QUARANTINE"}:
            raise ValueError(f"Unknown quality tier for {mutation}")
        records = grouped[mutation]
        if int(row["n_rows"]) != len(records):
            raise ValueError(f"Record count does not reconcile for {mutation}")
        expected_ids = set(filter(None, row["target_source_record_ids"].split("|")))
        if expected_ids != {r["record_id"] for r in records}:
            raise ValueError(f"Source IDs do not reconcile for {mutation}")
        contexts = defaultdict(list)
        for record in records:
            contexts[record["context_key"]].append(number(record["delta_log10_eff"]))
        verified = {r["context_key"] for r in records if r["context_status"] in {"PRESENT", "DOCUMENTED"}}
        if len(verified) != int(row["n_verified_contexts"]):
            raise ValueError(f"Contributor-group count does not reconcile for {mutation}")
        if records:
            observed = statistics.median(statistics.median(v) for v in contexts.values())
            if not math.isclose(observed, number(row["target_delta_log10_eff"]), abs_tol=1e-12):
                raise ValueError(f"Measured result does not reconstruct for {mutation}")
        if row["quality_tier"] == "QUARANTINE" and records:
            raise ValueError(f"Quarantined mutation has target records: {mutation}")
        kept = [r for r in records if not record_exclusion(r)]
        measured, context_values = summarize_records(kept)
        row_values = [number(r['delta_log10_eff']) for r in kept]
        threshold = policy['effect_threshold_log10']
        reversal = any(min(v) <= -threshold and max(v) >= threshold for v in (row_values, context_values) if v)
        high_spread = any(statistics.stdev(v) > policy['high_spread_threshold_log10'] for v in (row_values, context_values) if len(v)>1)
        repeated = len(records) >= policy["minimum_eligible_records"]
        enough = len(kept) >= policy['minimum_eligible_records']
        enough_groups = len(context_values) >= policy['minimum_contributor_groups']
        primary = enough and enough_groups and not reversal and not high_spread
        reasons = []
        if not enough:
            reasons.append(f"fewer_than_{policy['minimum_eligible_records']}_retained_records")
        if not enough_groups:
            reasons.append(f"fewer_than_{policy['minimum_contributor_groups']}_contributor_groups")
        if reversal:
            reasons.append('opposing_measured_effects')
        if high_spread:
            reasons.append('measurement_spread_above_existing_limit')
        direct = [r for r in kept if r['baseline_type']=='reference_wt']
        direct_value, _ = summarize_records(direct)
        ledger.append({
            "mutation_id": mutation, "position": int(row["position"]),
            "total_database_records": int(row.get('n_raw_rows', len(records))),
            "usable_records": len(records), "retained_records": len(kept),
            "contributor_groups": len(context_values),
            "quality_tier": row["quality_tier"],
            "all_repeated": repeated, "primary": primary,
            "highest_repeat": primary and len(kept)>=policy['highlight_minimum_records'],
            "retained_observed_log10": measured if measured is not None else '',
            "direct_wt_records": len(direct),
            "direct_wt_observed_log10": direct_value if direct_value is not None else '',
            "exclusion_reasons": "|".join(reasons),
            "existing_quality_flags": row["quality_reason_codes"],
            "source_record_ids": '|'.join(r['record_id'] for r in kept),
            "earlier_usable_record_ids": row["target_source_record_ids"],
            "earlier_excluded_record_ids": row.get('quarantined_record_ids', ''),
            "newly_excluded_usable_record_ids": '|'.join(r['record_id'] for r in records if record_exclusion(r)),
        })
    return sorted(ledger, key=lambda r: (-r['retained_records'], -r['total_database_records'], r["position"], r["mutation_id"]))


def ranks(values):
    ordered = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        for index in ordered[start:end]:
            result[index] = (start + end - 1) / 2 + 1
        start = end
    return result


def metrics(rows, threshold):
    if not rows:
        return {"n_mutations": 0, "spearman": None, "mean_absolute_error_log10": None,
                "no_change_mean_absolute_error_log10": None, "direction_correct": 0,
                "direction_total": 0, "beneficial_measured": 0, "beneficial_sign_correct": 0}
    observed = [r["observed_log10"] for r in rows]
    predicted = [r["predicted_log10"] for r in rows]
    correlation = None
    if len(rows) >= 2 and len(set(observed)) > 1 and len(set(predicted)) > 1:
        correlation = statistics.correlation(ranks(observed), ranks(predicted))
    directional = [r for r in rows if abs(r["observed_log10"]) >= threshold]
    beneficial = [r for r in rows if r["observed_log10"] >= threshold]
    return {
        "n_mutations": len(rows), "spearman": correlation,
        "mean_absolute_error_log10": statistics.mean(abs(a-b) for a,b in zip(observed,predicted)),
        "no_change_mean_absolute_error_log10": statistics.mean(abs(a) for a in observed),
        "direction_correct": sum(r["observed_log10"] * r["predicted_log10"] > 0 for r in directional),
        "direction_total": len(directional),
        "beneficial_measured": len(beneficial),
        "beneficial_sign_correct": sum(r["predicted_log10"] > 0 for r in beneficial),
    }


def compare(ledger, consensus, measurements, predictions, model, policy):
    """Require complete predictions; never drop a selected mutation silently."""
    chosen = [r for r in predictions if r["evaluation"] == "position_holdout" and r["model"] == model]
    pred = unique(chosen, "mutation_id")
    con = unique(consensus, "mutation_id")
    results = []
    for row in ledger:
        if not row["all_repeated"]:
            continue
        mutation = row["mutation_id"]
        if mutation not in pred:
            raise ValueError(f"Missing prediction for selected mutation {mutation}")
        original_observed = number(con[mutation]["target_delta_log10_eff"])
        if not math.isclose(number(pred[mutation]["observed_delta_log10_eff"]), original_observed, abs_tol=1e-12):
            raise ValueError(f"Prediction table uses a different measured target for {mutation}")
        predicted = number(pred[mutation]["prediction"])
        kept_ids = set(row['source_record_ids'].split('|'))
        values = [number(r['delta_log10_eff']) for r in measurements if r['record_id'] in kept_ids]
        observed = number(row['retained_observed_log10']) if values else original_observed
        results.append({**row, "observed_log10": observed, "predicted_log10": predicted,
                        "earlier_observed_log10": original_observed,
                        "observed_fold_vs_wt": 10**observed, "predicted_fold_vs_wt": 10**predicted,
                        "record_min_fold_vs_wt": 10**min(values) if values else '', "record_max_fold_vs_wt": 10**max(values) if values else ''})
    summaries = {key: metrics([r for r in results if r[key]], policy["effect_threshold_log10"])
                 for key in ("primary", "highest_repeat")}
    summaries['all_repeated'] = metrics([{**r, 'observed_log10':r['earlier_observed_log10']} for r in results], policy['effect_threshold_log10'])
    return results, summaries


def write_csv(path, rows):
    if not rows:
        raise ValueError("Cannot export an empty table")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def report_text(summary, comparison):
    title = "CatPred results" if summary['model'] == 'raw_catpred' else "model results"
    text = [f"# Repeated BglB measurements and {title}", "",
            f"Dataset: {summary['dataset_label']}.", "",
            "Each dot is a measurement included in this comparison. Different database entries do not establish independent biological replicates.",
            "Entries with the same recorded contributor are grouped together. These groups are not verified experimental batches.", "",
            "## Which measurements are included?", "",
            f"The input contains {summary['selection']['total_mutations']} substitutions. "
            f"{summary['selection']['all_repeated']} have at least four usable measurements. "
            f"{summary['selection']['primary']} meet the additional selection rules, and "
            f"{summary['selection']['highest_repeat']} of those have six or more measurements included. "
            "Mutations are ordered by the number of included measurements. Agreement with the model does not affect selection.", "",
            "An entry may be excluded because contributor information is missing, a suitable WT reference is unavailable, or a data issue was flagged in the earlier checks. The selection also limits opposing changes and large differences among measurements. An exclusion does not prove that a measurement is noise.", "",
            "The selection was revised after the broader results were known. The methods explain that change, including why a separate excluded entry no longer prevents using the other measurements for a mutation.", "",
            "## Main comparison", "",
            "| Mutation | All entries | Usable measurements | Measurements included | Contributor groups | Measured efficiency / WT | Predicted efficiency / WT | Range of included measurements |",
            "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in comparison:
        if row["primary"]:
            text.append(f"| {row['mutation_id']} | {row['total_database_records']} | {row['usable_records']} | {row['retained_records']} | {row['contributor_groups']} | "
                        f"{row['observed_fold_vs_wt']:.3g} | {row['predicted_fold_vs_wt']:.3g} | "
                        f"{row['record_min_fold_vs_wt']:.3g} to {row['record_max_fold_vs_wt']:.3g} |")
    text += ["", "WT means the wild-type enzyme. Both measured and predicted values show catalytic efficiency relative to WT. A value of 1 means no change; 0.1 means one tenth of WT.",
             "To calculate the combined measured value, express each mutant-to-WT ratio on a log10 scale, take the median within each contributor group and then the median across groups, and convert back to a ratio. This gives each contributor group equal weight. The range shows the lowest and highest included measurements; it is not a confidence interval.",
             "", "## Results for the three comparisons", "",
             "| Set | Mutations | Rank agreement | Average error (log10) | No-change average error (log10) | Correct direction |",
             "|---|---:|---:|---:|---:|---:|"]
    names = {"primary": "All mutations that meet the selection rules",
             "highest_repeat": "Mutations meeting the rules with at least six measurements included",
             "all_repeated": "All mutations with at least four usable measurements, using the original combined values"}
    for key, values in summary["metrics"].items():
        fmt = lambda value: "not available" if value is None else f"{value:.3f}"
        text.append(f"| {names[key]} | {values['n_mutations']} | {fmt(values['spearman'])} | "
                    f"{fmt(values['mean_absolute_error_log10'])} | {fmt(values['no_change_mean_absolute_error_log10'])} | "
                    f"{values['direction_correct']}/{values['direction_total']} |")
    text += ["", "Rank agreement is Spearman correlation: 1 means the same ranking and -1 means the reverse ranking. Small sets make this number unstable.",
             "Average error uses the log10 scale, so large fold differences do not dominate solely because of their units. The no-change comparison always predicts WT-like activity.",
             "Direction is scored only when the measured effect is at least 0.2 log10 from WT, using the prediction's sign. Tiny predicted effects can therefore count as correct direction.",
             "", "## Limits", "",
             "This comparison uses the BglB/pNPG data and model results identified in the input files. The model was not rerun for this report. Database dates and model settings are documented in the project methods.",
             "The broader dataset results were already known before this subset analysis. No untouched test set, new model fitting, or new experimental validation is claimed.",
             "Four entries alone do not establish reliability. The selection rules limit opposing changes and large differences among the included measurements. A mutation can still be included when its other measurements meet the rules.",
             "Selection remains retrospective, and the set is small. Results cannot establish performance across D2D, other enzymes, or other substrates.",
             "The earlier audit is required: this tool does not independently validate raw sequences, units, expression labels, or source data.", ""]
    return "\n".join(text)


def plot(comparison, measurements, output, model_label="Model result", highest_only=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    subset = 'highest_repeat' if highest_only else 'primary'
    rows = [r for r in comparison if r[subset]]
    if not rows:
        raise ValueError("No mutations in the primary set to plot")
    fig, ax = plt.subplots(figsize=(max(9, len(rows)*1.4), 5.8))
    for i, row in enumerate(rows):
        ids = set(row['source_record_ids'].split('|'))
        values = [10**number(r["delta_log10_eff"]) for r in measurements if r['record_id'] in ids]
        offsets = [(j-(len(values)-1)/2)*0.025 for j in range(len(values))]
        ax.scatter([i+v for v in offsets], values, s=32, color="#7b8993", alpha=0.8,
                   label="Included measurement" if i==0 else None)
        ax.scatter(i-0.18, row["observed_fold_vs_wt"], marker="D", color="#176c68", s=55,
                   label="Combined measured value" if i==0 else None)
        ax.scatter(i+0.18, row["predicted_fold_vs_wt"], marker="X", color="#c05b27", s=70,
                   label=model_label if i==0 else None)
    ax.axhline(1, color="#343e47", linestyle="--", linewidth=1, label="Wild-type enzyme (WT)")
    ax.set_yscale("log")
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.set_ylabel("Catalytic efficiency relative to WT (fold, log scale)")
    ax.set_xticks(range(len(rows)), [f"{r['mutation_id']}\n{r['retained_records']} included\n{r['total_database_records']} total entries" for r in rows])
    ax.set_title("BglB mutations with the most included measurements" if highest_only else "All BglB mutations that meet the selection rules", loc="left", pad=20)
    ax.grid(axis="y", alpha=0.15)
    ax.spines[["right", "top"]].set_visible(False)
    ax.legend(loc="lower left", fontsize=8)
    fig.text(0.08, 0.02, "Dots are included measurements; total entries also count missing or excluded measurements.\nCombined value: medians within and across contributor groups on the log10 scale, converted back to a ratio.\nSeparate entries do not establish independent experiments.", fontsize=9)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def run(args):
    sources = {k: Path(getattr(args, k)).resolve() for k in ("consensus", "measurements", "predictions", "policy")}
    policy = load_policy(sources["policy"])
    consensus = read_csv(sources["consensus"], ["mutation_id", "position", "n_rows", "n_verified_contexts", "quality_tier", "quality_reason_codes", "target_delta_log10_eff", "target_source_record_ids"])
    measurements = read_csv(sources["measurements"], ["record_id", "mutation_id", "context_key", "context_status", "baseline_type", "delta_log10_eff", "phase1_downweight_reason_codes"])
    ledger = select(consensus, measurements, policy)
    predictions = read_csv(sources["predictions"], ["mutation_id", "evaluation", "model", "prediction", "observed_delta_log10_eff"])
    comparison, scores = compare(ledger, consensus, measurements, predictions, args.model, policy)
    if not comparison:
        raise ValueError("No mutations have at least four usable measurements")
    output = Path(args.out).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output must be a new or empty directory; existing results are preserved")
    summary = {"dataset_label": args.dataset_label, "model": args.model, "policy": policy,
               "selection": {"total_mutations": len(ledger), **{k: sum(r[k] for r in ledger) for k in ("primary", "highest_repeat", "all_repeated")}},
               "metrics": scores,
               "software": {"version": "0.3.0", "python": sys.version.split()[0],
                            "analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
               "input_files": {k: {"name": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in sources.items()}}
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output/"selection.csv", ledger)
    write_csv(output/"comparison.csv", comparison)
    available = unique(measurements, 'record_id')
    by_mutation = unique(ledger, 'mutation_id')
    record_ledger = []
    for variant in consensus:
        for identifier in filter(None, variant.get('source_record_ids', variant['target_source_record_ids']).split('|')):
            record = available.get(identifier)
            reason = record_exclusion(record) if record else 'excluded_in_earlier_raw_data_audit'
            record_ledger.append({'mutation_id':variant['mutation_id'], 'record_id':identifier,
                                  'retained_for_summary':not bool(reason),
                                  'included_in_primary_comparison':not bool(reason) and by_mutation[variant['mutation_id']]['primary'],
                                  'row_exclusion_reason':reason})
    write_csv(output/'record_selection.csv', record_ledger)
    (output/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    (output/"report.md").write_text(report_text(summary, comparison))
    if args.plot:
        label = "CatPred result" if args.model == "raw_catpred" else f"{args.model} result"
        plot(comparison, measurements, output/"comparison.png", label, highest_only=bool(summary['selection']['highest_repeat']))
        plot(comparison, measurements, output/"all_selected.png", label)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("consensus", "measurements", "predictions", "policy", "out", "dataset-label"):
        parser.add_argument("--"+name, required=True)
    parser.add_argument("--model", default="raw_catpred")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    try:
        result = run(args)
    except (ValueError, KeyError, OSError, ImportError) as error:
        parser.exit(2, f"Cannot complete comparison: {error}\n")
    print(json.dumps(result["selection"], indent=2))
