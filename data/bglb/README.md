# Real analysis inputs

The experimental data were obtained from the public
[Design2Data database](https://d2d-cure.vercel.app/database/characterization_data/BglB).
The measurements belong to the work of D2D and its data contributors.
Mixin Zhao wrote the analysis code. CatPred was developed by its own authors.

These files are the numerical inputs for this particular BglB/pNPG comparison.
They preserve the original snapshot and audit decisions. They are not a new
download of the full current database.

| File | Contents |
|---|---|
| `audit.csv` | 1,519 records: measured values, reference values, sequence checks, and earlier inclusion or exclusion decisions |
| `measurements.csv` | 968 usable mutant records, their WT reference sources, and calculated changes relative to WT |
| `consensus.csv` | 756 individual mutations, original measured summaries, and all contributing or excluded record IDs |
| `predictions.csv` | CatPred results for 615 mutations with usable measured targets, including the values used to normalize against WT |
| `reference.json` | The 449-amino-acid WT sequence and the original coding-region decision |
| `manifest.json` | Source names and hashes, packaged-file hashes, dates, and known limits |

All 1,519 records remain represented, including 241 excluded records. Of
these exclusions, 224 are mutant records and 17 are WT records. The other
310 WT records provide reference values. Exclusion is not the same as zero
activity and does not prove that a measurement is noise.

## Read the fields

- `record_id` is the original public D2D record ID.
- `variant_norm` and `mutation_id` are uppercase mutation labels. The raw
  label `D150v` is represented as `D150V`; the underlying record ID is unchanged.
- `context_key` identifies the original contributor group. Names were replaced
  with sequential `group_0001`-style labels. These are not experiment or batch
  IDs. Unknown contributors remain unknown; each unknown eligible record has
  a separate `UNVERIFIED::record_id` grouping key.
- `kcat_value` is in 1/min and `km_value` is in mM, following the D2D export
  headers. Experimental efficiency is in 1/(mM min). Positive component values
  are divided to calculate efficiency. A positive reported efficiency is used
  only when the two component values are not both usable. Empty means missing.
- `baseline_type` is `reference_wt` for the WT values recorded with that row,
  `context_wt` for the median of usable WT records from that contributor group,
  or `global_wt` for the median of all usable WT records. Source record IDs
  are included. Global-WT comparisons remain available but are excluded from
  the selected set.
- `delta_log10_eff` is `log10(mutant efficiency / WT efficiency)`. Zero means
  WT-like efficiency; `10 ** delta_log10_eff` gives the displayed fold value.
- `target_delta_log10_eff` is the original mutation summary, calculated by
  taking medians within groups and then across groups on the log10 scale.
  The selected-set summary is recalculated after the additional row filters.
- `quality_tier` preserves the earlier audit labels. The revised selection
  evaluates individual rows and does not require the whole mutation to be CORE.
- `catpred_efficiency_ratio` and `catpred_wt_efficiency_ratio` are model-output
  quantities used for WT normalization, not experimental measurements. Their
  ratio gives the predicted fold change. `catpred_input_rows` counts prediction
  input rows, not independent experiments.

## Understand exclusions and flags

| Reason | Meaning |
|---|---|
| `MISSING_REQUIRED_FIELD` | Information required by the earlier audit was absent. The action depends on the missing field: some records were excluded; others remained usable with lower confidence. |
| `EXPRESSION_CONTRADICTION` | The expression information did not support using the record for kinetic regression. This was not converted to zero activity. |
| `EFFICIENCY_COMPONENT_CONFLICT` | Reported efficiency and efficiency calculated from kcat and KM disagreed. |
| `WT_BASELINE_GLOBAL_ONLY` | Only the overall WT reference was available. |
| `ID_JOIN_MISSING` | The record did not have a row matched by ID in the earlier CatPred input. Its expected sequence could still be generated from the mutation label. |

The original exclusion reasons and row flags are preserved. A flag does not
always mean exclusion: `phase1_downweight_reason_codes` records the particular
earlier decisions used by the additional row filter. The generated
`record_selection.csv` states both the earlier exclusion reason and any new
reason. The [methods](../../docs/methods.md) describe the remaining mutation
checks, including opposing effects and measurement spread.

## Provenance and reproduction

Run `bglb-reproduce --out results --plot` from the repository root after
installation. This checks every input hash, measured efficiency, WT baseline,
mutation sequence hash, source-record count, original measured summary, and
CatPred WT normalization before reporting the comparison.

The original audit decisions are supplied as inputs. The checks do not
recreate missing instrument records or verify experimentally that every
protein had the intended sequence. The exact code revision of the May 15,
2026 CatPred export is unknown; model inference is not rerun.

`tools/prepare_release_data.py` records how these files were extracted from
the earlier audit. It selects numerical and decision fields, replaces
contributor labels consistently, keeps record IDs, selects one CatPred result
per eligible mutation, and writes file hashes. It does not alter source files.
Its `--help` lists the source inputs needed to prepare another release.

The original download date is unknown. The source file's modification date is
February 3, 2026, and the audit was frozen in July 2026. A September 24, 2026
check found every source record through the public endpoint and matched all
nine kinetic fields in the original export. See the
[source-check record](../../docs/d2d_source_check.json) for changes to other
fields and current approval status.

The repository's MIT license covers its code and documentation. It does not
relicense D2D's experimental data or the CatPred model. This numerical extract
keeps source attribution and omits contributor names and institutions.
