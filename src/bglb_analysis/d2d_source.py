"""Check a local D2D export against the data used by the public database page."""

import argparse
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from urllib.request import urlopen

from .analysis import read_csv, unique


PAGE_URL = "https://d2d-cure.vercel.app/database/characterization_data/BglB"
DATA_URL = "https://d2d-cure.vercel.app/api/getCharacterizationData?enzyme=BglB"
FIELDS = {
    "Yield (mg/mL)": "yield_avg",
    "KM (mM)": "KM_avg", "KM SD": "KM_SD", "reference KM": "KM_ref",
    "kcat (1/min)": "kcat_avg", "kcat SD": "kcat_SD", "reference kcat": "kcat_ref",
    "kcat/KM (1/(mM min))": "kcat_over_KM", "kcat/KM SD": "kcat_over_KM_SD",
    "reference kcat/KM": "kcat_over_KM_ref",
    "T50 (degrees C)": "T50", "T50 SD": "T50_SD", "reference T50": "T50_ref",
    "Tm (degrees C)": "Tm", "Tm SD": "Tm_SD", "Rosetta score change": "Rosetta_score",
    "Institution": "institution", "Created by": "creator",
}
KINETIC_FIELDS = [k for k in FIELDS if "KM" in k or "kcat" in k]
AUDIT_KINETICS = {
    "km_value": "KM_avg", "kcat_value": "kcat_avg",
    "reported_efficiency": "kcat_over_KM", "reference_km_value": "KM_ref",
    "reference_kcat_value": "kcat_ref", "reported_reference_efficiency": "kcat_over_KM_ref",
}


def check_audit(local, public):
    """Compare the published numerical extract without needing creator names."""
    local_by_id = unique(local, "record_id")
    if not isinstance(public, list) or not public:
        raise ValueError("Public response must be a nonempty record list")
    current = unique([{**row, "id": str(row["id"])} for row in public], "id")
    missing, differences, uncurated, not_displayed = [], [], [], []
    field_counts = Counter()
    kinetic_mismatches = 0
    for identifier, row in local_by_id.items():
        if identifier not in current:
            missing.append(identifier)
            continue
        remote = current[identifier]
        for field in ("expressed", "curated", "submitted_for_curation"):
            if type(remote[field]) is not bool:
                raise ValueError(f"Public {field} must be boolean")
        changed = []
        for field, remote_field in AUDIT_KINETICS.items():
            a, b = row[field], remote[remote_field]
            # The original D2D CSV export represented numeric zero as blank.
            if not a and not b:
                continue
            if not a or b is None or not math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12):
                changed.append(field)
        kinetic_mismatches += bool(changed)
        variant = "WT" if remote["resid"] == "X" else f"{remote['resid']}{remote['resnum']}{remote['resmut']}".upper()
        if row["variant_norm"] != variant:
            changed.append("variant_norm")
        if row["expressed_yes"].lower() != str(remote["expressed"]).lower():
            changed.append("expressed_yes")
        if changed:
            differences.append({"record_id": identifier, "mutation": row["variant_norm"], "fields": changed})
            field_counts.update(changed)
        if not remote["curated"]:
            uncurated.append(identifier)
            if not remote["submitted_for_curation"]:
                not_displayed.append(identifier)
    return {
        "local_records": len(local), "public_endpoint_records": len(public),
        "local_ids_found": len(local) - len(missing), "missing_ids": missing,
        "local_records_matching_kinetic_columns": len(local) - len(missing) - kinetic_mismatches,
        "compared_columns": list(AUDIT_KINETICS) + ["variant_norm", "expressed_yes"],
        "differing_cells": sum(field_counts.values()), "differences_by_field": dict(field_counts),
        "differences": differences, "not_currently_curated_ids": uncurated,
        "not_currently_displayed_ids": not_displayed,
    }


def export_row(row):
    """Match the public page's CSV formatting, including its zero-to-blank rule.

    Induced? is a fixed export placeholder, not an induction-status check.
    This only formats values for comparison; it never edits the local table.
    """
    required = set(FIELDS.values()) | {"id", "resid", "resnum", "resmut", "expressed", "curated", "submitted_for_curation"}
    if required - row.keys():
        raise ValueError("Public record is missing required fields")
    for field in ("expressed", "curated", "submitted_for_curation"):
        if type(row[field]) is not bool:
            raise ValueError(f"Public {field} must be boolean")
    result = {field: str(row[key]) if row[key] else "" for field, key in FIELDS.items()}
    result.update({
        "ID #": str(row["id"]),
        "Variant": "WT" if row["resid"] == "X" else f"{row['resid']}{row['resnum']}{row['resmut']}",
        "Induced?": "data not transfered from old site",
        "Expressed?": "yes" if row["expressed"] else "no",
        "Curated?": "yes" if row["curated"] else "no",
        "Created by": row["creator"] or "unknown",
        "Yield (mg/mL)": str(row["yield_avg"]) if row["yield_avg"] is not None else ("not reported" if row["expressed"] else ""),
    })
    return result


def equivalent(left, right):
    if left == right:
        return True
    try:
        a, b = Decimal(left), Decimal(right)
        return a.is_finite() and b.is_finite() and a == b
    except InvalidOperation:
        return False


def check_records(local, public, comparison=()):
    local_by_id = unique(local, "ID #")
    if not isinstance(public, list) or not public:
        raise ValueError("Public response must be a nonempty record list")
    public_by_id = unique([{**r, "id": str(r['id'])} for r in public], "id")
    # New, unrelated entries can be incomplete. Only validate/export the
    # records actually used by this local dataset.
    formatted = {key: export_row(row) for key, row in public_by_id.items() if key in local_by_id}
    columns = set(FIELDS) | {"ID #", "Variant", "Induced?", "Expressed?", "Curated?"}
    if not local or any(set(row) != columns for row in local):
        raise ValueError("Local export must have the documented 23 columns")
    missing, differences, uncurated, not_displayed = [], [], [], []
    field_counts = Counter()
    for identifier, row in local_by_id.items():
        if identifier not in public_by_id:
            missing.append(identifier)
            continue
        current = public_by_id[identifier]
        changed = sorted(k for k in columns if not equivalent(row[k], formatted[identifier][k]))
        if changed:
            differences.append({"record_id": identifier, "mutation": row["Variant"], "fields": changed})
            field_counts.update(changed)
        if not current["curated"]:
            uncurated.append(identifier)
            if not current["submitted_for_curation"]:
                not_displayed.append(identifier)
    differs = {r['record_id'] for r in differences}
    kinetic_differs = {r['record_id'] for r in differences if set(r['fields']) & set(KINETIC_FIELDS)}
    sets = {}
    for name in ("primary", "highest_repeat", "all_repeated"):
        chosen = [r for r in comparison if r[name] == "True"]
        if not chosen:
            continue
        mutations = {r['mutation_id'] for r in chosen}
        ids = {identifier for identifier, row in local_by_id.items() if row['Variant'].strip().upper() in mutations}
        retained = {i for row in chosen for i in row['source_record_ids'].split('|') if i}
        if not retained <= ids:
            raise ValueError("Comparison refers to records outside the local mutations")
        sets[name] = {
            "mutations": len(mutations), "total_local_records": len(ids),
            "retained_records": len(retained),
            "missing_ids": sorted(ids & set(missing)),
            "differing_ids": sorted(ids & differs),
            "not_currently_curated_ids": sorted(ids & set(uncurated)),
            "not_currently_displayed_ids": sorted(ids & set(not_displayed)),
        }
    return {
        "local_records": len(local), "public_endpoint_records": len(public),
        "public_curated_records": sum(r['curated'] for r in public),
        "public_displayable_records_including_pending": sum(r['curated'] or r['submitted_for_curation'] for r in public),
        "local_ids_found": len(local) - len(missing), "missing_ids": missing,
        "local_records_matching_all_columns": len(local) - len(missing) - len(differs),
        "local_records_matching_kinetic_columns": len(local) - len(missing) - len(kinetic_differs),
        "compared_columns": sorted(columns), "kinetic_columns": KINETIC_FIELDS,
        "differing_cells": sum(field_counts.values()), "differences_by_field": dict(sorted(field_counts.items())),
        "differences": differences, "not_currently_curated_ids": uncurated,
        "not_currently_displayed_ids": not_displayed, "comparison_sets": sets,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--database", type=Path, help="Original 23-column D2D export")
    source.add_argument("--audit", type=Path, help="Bundled numerical audit.csv")
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("Output already exists; choose a new path")
    if args.audit and args.comparison:
        parser.error("--comparison is only supported with --database")
    input_path = args.audit or args.database
    local = read_csv(input_path, ["record_id", "variant_norm", "expressed_yes", *AUDIT_KINETICS] if args.audit else ["ID #", "Variant"])
    comparison = read_csv(args.comparison, ["mutation_id", "primary", "highest_repeat", "all_repeated", "source_record_ids"]) if args.comparison else []
    # No login, cookies, API key, or other credentials are sent.
    with urlopen(DATA_URL, timeout=60) as response:
        payload, status = response.read(), response.status
    summary = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "database_page": PAGE_URL, "public_endpoint": DATA_URL,
        "authentication": "none", "http_status": status,
        "public_response_sha256": hashlib.sha256(payload).hexdigest(),
        "public_response_bytes": len(payload),
        "local_database_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "comparison_sha256": hashlib.sha256(args.comparison.read_bytes()).hexdigest() if args.comparison else None,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "method": ("Join by record ID; compare six numerical fields with relative and absolute tolerance 1e-12, mutation labels, and expression flags. Preserve the original export's zero-to-blank convention. No data values are changed." if args.audit else "Join by record ID; compare all 23 export columns. Numbers must be exactly equal after decimal parsing. Follow the public page's CSV blank and WT-label formatting. No data values are changed."),
        "scope": "Endpoint availability is recorded separately from current table display and curation. This check does not establish the original download date, permission to redistribute the database, or experimental independence.",
        **(check_audit(local, json.loads(payload)) if args.audit else check_records(local, json.loads(payload), comparison)),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({k: summary[k] for k in ('local_records', 'local_ids_found', 'local_records_matching_kinetic_columns', 'differing_cells')}, indent=2))


if __name__ == "__main__":
    main()
