# BglB mutation analysis

[![Tests](https://github.com/MixinZh/bglb-mutation-analysis/actions/workflows/tests.yml/badge.svg)](https://github.com/MixinZh/bglb-mutation-analysis/actions/workflows/tests.yml)

**Do CatPred results match repeatedly measured changes in BglB catalytic efficiency?**

BglB is a beta-glucosidase, an enzyme that breaks certain sugar bonds. This
project compares CatPred results with repeated measurements of individual
BglB mutations. The comparison asks whether CatPred captures changes in
catalytic efficiency, a measure of how effectively the enzyme processes its
substrate, relative to the wild-type enzyme (WT).

**Data source:** The experimental data were obtained from the public
[Design2Data (D2D) database](https://d2d-cure.vercel.app/database/characterization_data/BglB).
Credit for the experimental measurements belongs to D2D and its data contributors.

## What did the comparison find?

**CatPred missed some repeatedly observed BglB changes, while matching others.**
For example, all 11 N220F measurements included here show higher catalytic
efficiency than WT, whereas the CatPred result is below WT. H315N provides
an agreement case: both its combined measured value and its CatPred result
are close to WT.

The figure highlights four mutations with at least six measurements included
in this comparison. Among the mutations that meet the selection rules, these
have the most included measurements. Agreement with CatPred does not affect
which mutations are shown.

![Repeated BglB measurements compared with CatPred results](docs/figures/comparison.png)

Both the measured values and CatPred results below show catalytic efficiency
relative to WT. **WT is 1:** a value above 1 means higher efficiency, and a
value below 1 means lower efficiency.

| Mutation | All database entries | Entries with usable measurements | Measurements included | Measured efficiency / WT | CatPred result / WT |
|---|---:|---:|---:|---:|---:|
| N220F | 12 | 11 | 11 | 3.36 | 0.58 |
| H315N | 15 | 8 | 7 | 1.04 | 1.03 |
| R246K | 13 | 7 | 7 | 1.82 | 0.95 |
| N404M | 8 | 7 | 6 | 0.0075 | 0.62 |

The **combined measured value** is calculated in three steps. First, express
each measurement's mutant-to-WT ratio on a log10 scale. Then take the median
within each group of entries with the same recorded contributor, followed by
the median across those groups. Finally, convert that value back to a ratio
relative to WT. This gives each contributor group equal weight. It is not an
arithmetic average of all the measurements.

All database entries are counted in the first column, including entries
without usable kinetic measurements. Some entries with usable measurements
are also left out because they lack a suitable WT reference, contributor
information, or meet another stated exclusion rule. The following section
explains the selection.

For **N220F**, the 11 included measurements range from 1.21 to 12.68 times WT.
Their size varies, but every value is above WT. The combined measured value
is 3.36 times WT, compared with a CatPred result of 0.58 times WT. Repeating
the comparison using only the ten entries with their own recorded WT
reference values still leaves every measured value above WT.

For **R246K**, the combined measured value is above WT, but the individual
measurements range from 0.68 to 4.00 times WT. The direction is less consistent
than for N220F. For **N404M**, all six included measurements are below 0.012
times WT, while the CatPred result is 0.62 times WT.

## How were measurements chosen?

1. Start with entries that passed the earlier checks for mutation identity,
   expression, kinetic measurements, and units.
2. Include measurements with a recorded WT reference or WT measurements from
   the same contributor group, known contributor information, and no data
   issue that the earlier checks marked as reducing confidence in that measurement.
3. Require at least four included measurements from at least four recorded
   contributor groups. Apply the existing limits on measurements that show
   opposing changes or a large spread in values.
4. Order the mutations by the number of included measurements.

**Eleven mutations meet these rules.** In addition to the four above, they
are N270W, W325Q, C38S, E406Q, H223A, W325M, and E340A.
[View the comparison for all 11 mutations](docs/figures/all_selected.png).

A broader comparison includes all **33 mutations with at least four usable
measurements**, before the additional selection above. It keeps the original
combined measured values and gives more mixed results. Both comparisons are
reported by the code. The [methods](docs/methods.md#revised-selection) explain
the exact limits, calculation rules, and earlier change to the selection.
Every entry's inclusion or exclusion is recorded in the output files.

## What are the limits?

Different database entries and contributor groups do **not** establish
independent biological replicates. N220F entries 616 and 693 have identical
kcat and KM values but different WT-reference information. Whether they
come from the same experiment is unresolved. The comparison using only
entries with their own WT reference excludes entry 616 and still finds
higher efficiency than WT in every remaining N220F measurement.

An excluded entry is not automatically noise. Missing information,
expression failure, and conflicting measurements are different issues.
A mutation can still be included when its other measurements meet the rules.

The selection was revised after the broader results were known. This is a
retrospective comparison, not an untouched test set. It covers the BglB
mutations examined here using pNPG. It does not establish performance across
all of D2D, other enzymes, or other substrates.

The CatPred results come from a May 2026 export whose exact code revision is
unknown. This project reproduces the comparison from those results; it does
not rerun CatPred or test its current release. The D2D inputs are also an
earlier snapshot, rather than the full current database.

## Reproduce the real results

Python 3.10 or later is required. From a fresh clone:

```bash
git clone https://github.com/MixinZh/bglb-mutation-analysis.git
cd bglb-mutation-analysis
python -m venv .venv
source .venv/bin/activate
python -m pip install '.[plots]'
bglb-reproduce --out results --plot
python -m unittest discover -s tests -v
```

On Windows, activate the environment with `.venv\Scripts\activate`.
The real numerical inputs are included in [data/bglb](data/bglb). No lab
folder, account, GPU, or model download is needed. The calculation itself
uses the Python standard library; Matplotlib is only needed for figures.
For tables alone, install `.` and omit `--plot`.

The command checks the input file hashes, recalculates the measured changes
from the kinetic values and WT references, checks mutation sequences against
WT, and verifies the normalization of the CatPred results. It then recreates
the selection, tables, figures, and metrics. Tests compare the results with
the earlier verified analysis. Existing nonempty output folders are preserved;
use another `--out` folder for a new run.

Outputs include `selection.csv`, `record_selection.csv`, `comparison.csv`,
`summary.json`, and `report.md`. With plotting enabled, `comparison.png`
shows the four mutations with the most included measurements, and
`all_selected.png` shows all 11 mutations that meet the selection rules.

## Use the comparison code with other inputs

An invented example shows the required table format:

```bash
bglb-compare --consensus examples/synthetic/consensus.csv --measurements examples/synthetic/measurements.csv --predictions examples/synthetic/predictions.csv --policy selection_policy.json --model demo_model --dataset-label SYNTHETIC-DEMO --out demo_results
```

This example provides no evidence about CatPred. The generic command accepts
other audited tables in the same format. It checks record counts, grouping,
and combined measured values before comparing them with model results. The
[input guide](docs/methods.md#required-inputs) describes the fields.

## Data and software included

The repository includes numerical measurements and earlier data-check
results for all 1,519 source entries, 968 usable mutant measurements, a table
covering 756 mutations, and CatPred results for the 615 mutations
with usable measured values for comparison. Contributor names are replaced
with group IDs. Names, institutions, lab notebooks, and instrument files are
not included. Each entry keeps its D2D ID so readers can check it against the
public database.

On September 24, 2026, all 1,519 source entries were found through D2D's public
data endpoint, and all kinetic values matched. The
[methods](docs/methods.md#public-source-check) and
[source-check results](docs/d2d_source_check.json) describe changes to other
fields and approval status. The [data guide](data/bglb/README.md) explains
the input files and their dates.

The model results were generated with
[CatPred](https://github.com/maranasgroup/CatPred), developed by its authors.
This is not an official D2D or CatPred repository. The code reads previously
checked inputs; it does not train models or recreate the original raw-data
checks.

## Cite or contribute

Use [CITATION.cff](CITATION.cff) to cite this analysis. Please also credit
[D2D](https://d2d.ucdavis.edu/d2d-database-app) for the measurements and
[CatPred](https://github.com/maranasgroup/CatPred) for the prediction model.
Corrections and improvements are welcome through GitHub issues or pull
requests. Include a small example, the expected behavior, and the test result.

## License

The code is copyright 2026 Mixin Zhao and is released under the
[MIT License](LICENSE). This license covers this repository's code and
documentation. D2D data and CatPred remain subject to their own terms.
