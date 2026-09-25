# Methods and inputs

## Authorship and data source

Mixin Zhao is the sole author and code contributor for this repository.
AI tools were used during code development and review.

The experimental data were obtained from the public
[Design2Data (D2D) database](https://d2d-cure.vercel.app/database/characterization_data/BglB).
Credit belongs to D2D and the people who generated the measurements.
[CatPred](https://github.com/maranasgroup/CatPred) is the model used to generate
the predictions. D2D and CatPred are external sources, not code contributors
to this repository.

## The unit is one named mutation

N220F means one amino-acid substitution. Records for other mutations at
position N220 are not counted as repeats of N220F. Different database records
do not necessarily come from separate experiments. The source file has
separate `Institution` and `Created by` columns. This analysis uses
`Created by`, the submitter username, to group records. A username does not
identify an institution or prove who performed the experiment.

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

The counts describe three different things:

- **D2D records:** All records for that mutation in the saved database, including records without usable measurements.
- **Usable measurements:** Measurements that passed the earlier data checks.
- **Included measurements:** Usable measurements that also pass the additional checks below.

Keep a measurement only if it has a known submitter username, a WT reference
saved with that record or calculated from WT measurements under the same
username, and no earlier flag in `phase1_downweight_reason_codes`.
Measurements that use only the database-wide WT reference are left out.
Their usable values remain in the broader comparison.

For example, N404M has seven usable measurements. Record 3330 uses the
database-wide WT reference, so it is excluded from the selected comparison.
The other six have their own recorded WT references and are included.

A mutation must have at least four included measurements from at least four
different submitter usernames. Exclude a mutation if its measurements strongly
conflict: at least one effect is at or below -0.2 log10 and another is at or
above +0.2 log10, about 0.63 and 1.58 times WT. Also exclude it if the sample
standard deviation exceeds 0.60 log10. Apply both checks to the individual
measurement effects and to the median effects for each username. These rules
do not require every measurement to be on the same side of WT. Exclusion
does not prove that a measurement is wrong.

To calculate the combined measured value, divide each mutant's efficiency
by its experimental WT reference and take log10. Take the median for each
username, then the median of those values. Convert back with `10 ** effect`.
This gives each username one value in the final median.

For the four mutations in the main figure, every included measurement has
a different username: N220F has 11 measurements from 11 usernames, H315N
7 from 7, R246K 7 from 7, and N404M 6 from 6. Combining measurements within
a username therefore has no effect on these four results. Their combined
value is simply the median of the included log10 ratios, converted back
to a ratio. This is not an arithmetic average of the ratios.

Order mutations that meet these rules by the number of included measurements,
then the total number of database records, then position and mutation name.
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
results were already known. The selection was changed after seeing the results, so this is not a test
on previously unseen results or a plan fixed before the comparison.

## Independence and the comparison using recorded WT references

N220F records 616 and 693 have identical kcat and KM values but different
WT-reference information. Their relationship is unresolved. Record 616 uses
WT measurements from the same submitter username; the other ten included
records have their own recorded WT reference values. Repeating the
comparison with only those ten records excludes 616 and still places every
measured result above WT. For every mutation, the output reports the number
of records with their own WT reference values and their combined measured
value.

R427L has four usable records but two submitter usernames, including potentially
reused kinetic records. It remains outside the selected set. Neither records
nor usernames prove that experiments were performed separately. That requires
records showing how each experiment was performed.

## Prediction comparison

Use `raw_catpred` rows from the earlier baseline table and its
`position_holdout` records to obtain one CatPred result per mutation.
CatPred itself was not fitted on those folds. Before comparing to the combined
measured value after the additional selection, verify that the prediction
table refers to the original measured value. Missing predictions stop the calculation.

The CatPred results come from the May 15, 2026 prediction export, which used
`log10kcat_max` and `log10km_mean`. Its exact code revision is not established
by the later adapter's pinned checkpoints. CatPred was not rerun for this
report, and the results do not establish performance of its current default
setup. The later adapter and E154 experiments are outside this comparison.

The comparison checks three things:

- **Order:** Does CatPred rank mutations by efficiency in the same order as the measurements? Spearman correlation measures this agreement.
- **Error:** How far are predictions from measurements? Mean absolute error is the average distance between their log10 ratios.
- **Increase or decrease:** Are the predicted and measured ratios both above 1 or both below 1? Only measured effects at least 0.2 log10 away from WT count, about 0.63 times WT or lower, or 1.58 times WT or higher. This cutoff applies to measurements, not predictions, and is not a test of statistical significance. Even a small predicted change can count as the correct direction.

The predicted ratio is predicted mutant efficiency divided by predicted WT
efficiency. The measured ratio uses the experimental WT reference. Each
mutation counts equally in these checks, regardless of its measurement count.
The code also compares CatPred with always predicting 1 times WT, meaning no
change. No statistical significance or confidence-interval claims are made.

The broad 33-mutation comparison keeps its original combined measured values.
Results for the 11 selected mutations use the combined values recalculated
from the measurements included after the additional selection.
Of those 33 mutations, 22 fail one or more additional checks: too few
measurements or usernames remain, measurements strongly conflict, or their
spread is too large. Disagreement with CatPred is not an exclusion rule.
Both comparisons are reported so readers can see the effect of selection.

## Required inputs

| Table | Required fields |
|---|---|
| Consensus | mutation_id, position, n_rows, n_verified_contexts, quality_tier, quality_reason_codes, target_delta_log10_eff, target_source_record_ids |
| Measurements | record_id, mutation_id, context_key, context_status, baseline_type, delta_log10_eff, phase1_downweight_reason_codes |
| Predictions | mutation_id, evaluation, model, prediction, observed_delta_log10_eff |

The full lab consensus also supplies `n_raw_rows`, `source_record_ids`, and
`quarantined_record_ids`, preserving total counts and earlier exclusions.
When absent, the minimal synthetic format accounts only for supplied records.
Consensus IDs are pipe-separated. Submitter username status is PRESENT or
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
for this calculation. No author-specific path is required. Submitter
usernames were replaced with consistent group IDs without changing grouping.

Before calculating the comparison, the command checks file hashes and
rebuilds the efficiencies, reference choices, and measured changes. It also
checks that source record counts reconcile, mutation sequences match the WT
reference, and CatPred values reproduce their reported WT-normalized results.
The original measured targets must reproduce from the grouped measurements.

`selection.csv` covers every mutation. `record_selection.csv` accounts for
all 1,192 mutant records, with specific earlier exclusion reasons.
`comparison.csv` retains all 33 repeated mutations with original and revised
summaries. `summary.json` records rules, metrics, validation counts, and input
and code hashes. Figures use the IDs of the included records, so excluded measurements
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
- All 86 records for the 11 selected mutations matched in every column and
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
submitter usernames or raw measurement values. That original check compared
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
