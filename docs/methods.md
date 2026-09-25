# Methods and inputs

## Authorship and data source

Mixin Zhao is the sole author and code contributor for this repository.
AI tools were used during code development and review.

The experimental data were obtained from the public
[Design2Data (D2D) database](https://d2d-cure.vercel.app/database/characterization_data/BglB).
D2D and its data contributors receive credit for those measurements.
[CatPred](https://github.com/maranasgroup/CatPred) is the model used to generate
the predictions. D2D and CatPred are external sources, not code contributors
to this repository.

## The unit is one named mutation

N220F means one amino-acid substitution. Records for other mutations at
position N220 are not counted as repeats of N220F. Different database records
and contributors do not establish independent biological experiments. The
earlier audit grouped records using the database's `Created by` field.

The raw snapshot has 1,519 records and 756 unique single substitutions. Its
recorded modification time is February 3, 2026; it was frozen in July 2026.
These dates do not establish when every experiment occurred. The public
database was checked on September 24, 2026. The original analysis inputs were
preserved; the check did not add new records or change existing values.

## Revised selection

Start with the measurements that passed the earlier data checks, included
in `data/bglb/measurements.csv`. Those checks covered mutation identity,
expression, usable positive kinetics, units, and WT comparisons. The bundled
inputs preserve the decisions. The reproduction command checks expected
sequences against WT and recalculates
efficiencies and WT comparisons. It cannot establish experimental sequence
identity or repeat checks requiring original instrument records.

Include a measurement only when it has WT reference values recorded with
that entry or calculated from WT measurements from the same contributor
group. Contributor information must be present, and the earlier checks must
not have marked a data issue as reducing confidence in that measurement.
In the input table, those issues appear in `phase1_downweight_reason_codes`.
Entries using only the overall WT reference, entries with unknown
contributors, and entries with those data issues are left out of this
comparison. Their original usable values remain in the broader comparison.

Calculate the combined measured value from the included measurements: take
the median within each contributor group and then the median across groups
on the `log10(mutant efficiency / WT efficiency)` scale. Convert back to the
displayed ratio with `10 ** effect`.

Require at least four included measurements and four contributor groups.
Apply the unchanged consistency checks to both the individual measurements
and the contributor-group medians. Leave a mutation out of the selected set
if either has values at or below -0.2 and at or above +0.2 log10, or sample
standard deviation above 0.60 log10. An exclusion does not prove noise;
passing these checks does not prove correctness.

Order mutations that meet these rules by the number of included measurements,
then the total number of database entries, then position and mutation name.
The main figure shows mutations with at least six included measurements.
This threshold controls presentation only: all mutations meeting the rules
with at least four included measurements are evaluated and shown in the full
figure. Agreement with CatPred does not affect selection or ordering.

## Correction to the earlier draft

The first draft required the old whole-mutation CORE label. N220F was excluded
because one of its 12 raw records was marked not expressed and had no usable
kinetics, even though its remaining 11 measurements were eligible. The revised
rule evaluates contributing records and does not allow a separate excluded
record to veto usable evidence. The same rule applies to every mutation.

Original audit labels are retained in the bundled inputs. Earlier model
results were already known. This is a retrospective correction, not an
untouched test set or preregistration.

## Independence and the comparison using recorded WT references

N220F records 616 and 693 have identical kcat and KM values but different
WT-reference information. Their relationship is unresolved. Record 616 uses
WT measurements from the same contributor group; the other ten included
entries have their own recorded WT reference values. Repeating the
comparison with only those ten entries excludes 616 and still places every
measured result above WT. For every mutation, the output reports the number
of entries with their own WT reference values and their combined measured
value.

R427L has four usable records but two contributor groups, including potentially
reused kinetic records. It remains outside the selected set. Neither records
nor contributors should be called independent biological replicates without
experiment-level evidence.

## Prediction comparison

Use `raw_catpred` rows from the earlier baseline table and its
`position_holdout` entries to obtain one CatPred result per mutation.
CatPred itself was not fitted on those folds. Before comparing to the combined
measured value after the additional selection, verify that the prediction
table refers to the original measured value. Missing predictions stop the calculation.

The CatPred results come from the May 15, 2026 prediction export, which used
`log10kcat_max` and `log10km_mean`. Its exact code revision is not established
by the later adapter's pinned checkpoints. CatPred was not rerun for this
report, and the results do not establish performance of its current default
setup. The later adapter and E154 experiments are outside this comparison.

Each mutation has equal weight. Metrics include rank agreement (Spearman
correlation), average absolute error on the log10 scale, and predicted direction
for measured effects at least 0.2 log10 away from WT. A no-change baseline
always predicts WT-like activity. No significance or confidence-interval
claims are made from this small set.

The broad 33-mutation comparison keeps its original combined measured values.
Results for the 11 selected mutations use the combined values recalculated
from the measurements included after the additional selection.
This distinction prevents the revised filtering from silently changing the
broader reference result.

## Required inputs

| Table | Required fields |
|---|---|
| Consensus | mutation_id, position, n_rows, n_verified_contexts, quality_tier, quality_reason_codes, target_delta_log10_eff, target_source_record_ids |
| Measurements | record_id, mutation_id, context_key, context_status, baseline_type, delta_log10_eff, phase1_downweight_reason_codes |
| Predictions | mutation_id, evaluation, model, prediction, observed_delta_log10_eff |

The full lab consensus also supplies `n_raw_rows`, `source_record_ids`, and
`quarantined_record_ids`, preserving total counts and earlier exclusions.
When absent, the minimal synthetic format accounts only for supplied records.
Consensus IDs are pipe-separated. Known contributor status is PRESENT or
DOCUMENTED; unknown status is PLACEHOLDER_UNKNOWN. Counts, source IDs, and
original summaries must reconcile before a revised calculation is allowed.

## Reproduce the real analysis

From the repository root after `python -m pip install '.[plots]'`:

```bash
bglb-reproduce --out results --plot
python -m unittest discover -s tests -v
```

The [bundled data](../data/bglb/README.md) contain the numerical measurements,
WT reference sources, original audit decisions, and CatPred results needed
for this calculation. No author-specific path is required. Contributor
names were replaced with consistent group IDs without changing grouping.

Before calculating the comparison, the command checks file hashes and
rebuilds the efficiencies, reference choices, and measured changes. It also
checks that source record counts reconcile, mutation sequences match the WT
reference, and CatPred values reproduce their reported WT-normalized results.
The original measured targets must reproduce from the grouped measurements.

`selection.csv` covers every mutation. `record_selection.csv` accounts for
all 1,192 mutant records, with specific earlier exclusion reasons.
`comparison.csv` retains all 33 repeated mutations with original and revised
summaries. `summary.json` records rules, metrics, validation counts, and input
and code hashes. Figures use the IDs of the included entries, so excluded measurements
cannot appear silently as included values. `report.md` provides a readable
summary. All outputs are written into the chosen output directory.

The tests compare the new outputs with the verified results made before
packaging. They also check that changed measurements, incorrect WT references,
missing predictions, and unexplained exclusions stop the calculation.
GitHub Actions runs the tests and both real and synthetic examples on
Python 3.10 and 3.13. Figure appearance can vary slightly with plotting-library
versions; numeric comparisons use a tolerance of 1e-12 where needed.

The original raw database hash is
`fdfa0a47c0ed9154532634e81bf1e9c3fd6302a24d756f07eac76bab44bffcee`.
The input manifest records hashes of the original audit and prediction tables
as well as hashes of their published numerical extracts. The original source
files were not changed.

## Public source check

D2D's [official database page](https://d2d.ucdavis.edu/d2d-database-app)
describes the database as a resource for testing enzyme prediction tools.
The project analyzes those measurements, with credit to D2D and the people
who generated them.

On September 24, 2026, a request without a login or credentials retrieved the
data used by the public database page. Records were joined by ID and all 23
columns in the local export were compared. Numeric values were compared
exactly after decimal parsing. The comparison follows the page's CSV
formatting: WT labels, blank fields, and the fixed `Induced?` placeholder.
Matching that placeholder does not verify induction status.
Subset membership uses the earlier audit's uppercase mutation labels; the
raw field comparison preserves the source spelling, including `D150v`.

- All 1,519 local record IDs were present, and all nine kinetic columns matched.
- 1,501 records matched in every column. The other 18 had 25 changed fields:
  14 approval flags and 11 temperature-related values.
- All 86 entries for the 11 selected mutations matched in every column and
  remain approved. This includes all 63 measurements in the selected comparison.
- In the broader 33-mutation set, F243H records 2575 and 1221 are no longer
  approved. Record 2575 is shown when pending records are included. Record
  1221 is returned by the data endpoint but is not currently displayed in
  either the approved or pending table. Keep that distinction visible;
  endpoint availability is not the same as current table inclusion.

The broader comparison retains the original snapshot and its labels. No
record was silently removed or updated. The current database contains more
records than the original snapshot, so this check does not make the analysis
a complete evaluation of today's database. The original download date remains
unknown; July 2026 is the freeze date.

The [machine-readable check](d2d_source_check.json) includes URLs, input and
response hashes, field counts, and record IDs for differences. It contains no
contributor names or raw measurement values. That original check compared
all 23 columns of the author's source export, including metadata omitted
from the published numerical extract.

To check the bundled numerical fields against the current public endpoint:

```bash
python -m bglb_analysis.d2d_source --audit data/bglb/audit.csv --out results/current_source_check.json
```

This command sends no login, cookies, or credentials. It checks the six
kinetic fields retained in the extract, mutation labels, and expression flags,
and reports current approval and display status. Small float-formatting
differences are allowed at a tolerance of 1e-12. It does not update the inputs.

For a full 23-column export obtained separately from D2D, the same checker
accepts `--database path/to/BglB_export.csv` instead of `--audit` and compares
numbers exactly after decimal parsing. Live results may change as D2D is
updated. Neither check establishes experimental independence.

## License and limits

The repository's code and documentation are covered by the
[MIT License](../LICENSE), copyright 2026 Mixin Zhao. The numerical inputs
are attributed to the public D2D database. The code license does not relicense
D2D data or CatPred.

This release supports reproducing and reviewing the stated comparison.
It does not resolve the original CatPred revision, confirm independent
experiments, evaluate the entire current database, or establish a generally
valid mutation predictor. Those limits do not prevent readers from checking
the calculations and the bounded findings reported here.
