# BglB mutation analysis

[![Tests](https://github.com/MixinZh/bglb-mutation-analysis/actions/workflows/tests.yml/badge.svg)](https://github.com/MixinZh/bglb-mutation-analysis/actions/workflows/tests.yml)

**Do CatPred results match enzyme changes that were measured repeatedly?**

BglB is a beta-glucosidase, an enzyme that breaks certain sugar bonds. This
project compares repeated measurements of individual BglB mutations with
CatPred results. It examines catalytic efficiency, a measure of how
effectively the enzyme processes the substrate, relative to the unchanged
enzyme (WT). 


**Data source:** The experimental data were obtained from the public
[Design2Data (D2D) database](https://d2d-cure.vercel.app/database/characterization_data/BglB).
Credit for the experimental measurements belongs to D2D and its data contributors.

## Start with the most repeated mutations

The lead figure shows qualifying mutations with six or more retained records.
All qualifying mutations with four or more retained records remain in the
full analysis.

![Repeated BglB measurements compared with CatPred results](docs/figures/comparison.png)

| Mutation | Total database records | Usable kinetic records | Retained for comparison | Measured fold vs WT | CatPred result |
|---|---:|---:|---:|---:|---:|
| N220F | 12 | 11 | 11 | 3.36 | 0.58 |
| H315N | 15 | 8 | 7 | 1.04 | 1.03 |
| R246K | 13 | 7 | 7 | 1.82 | 0.95 |
| N404M | 8 | 7 | 6 | 0.0075 | 0.62 |

WT is 1. Above 1 means higher catalytic efficiency; below 1 means lower
efficiency. Records without usable kinetics are not counted as repeats.
Some usable records are excluded from this comparison because they lack a
suitable WT reference or other required information.

**N220F provides a clear repeated increase.** All 11 retained values are above
WT, ranging from 1.21 to 12.68 times WT. Their summary is 3.36 times WT, while
the CatPred result is 0.58 times WT. The measured size varies, but the
direction is consistent. Restricting the check to the ten records with their
own WT reference still leaves every value above WT.

**H315N is an agreement case.** Its measured summary and prediction are both
close to WT. It stays in the figure because selection is based on measurement
support, not on how poorly CatPred performs.

R246K has a summary above WT, but its retained values range from 0.68 to 4.00
times WT. Its individual results are less uniform than N220F's. All six
retained N404M measurements are below 0.012 times WT; the prediction is 0.62.

In this comparison, **CatPred missed some repeatedly observed BglB changes,
while matching others.** The result applies to the mutations and measurements
examined here.

## How records are selected

1. Start with kinetic records that passed the earlier raw-data audit.
2. Retain records with a suitable WT reference, a known contributor, and no
   remaining row-level quality flag from that audit.
3. Require four retained records and four recorded contributor groups, and
   apply the existing checks for opposing effects and large measurement spread.
4. Rank qualifying mutations by retained-record count.

A separately excluded record does not erase usable evidence from the same
mutation. This corrects the first draft, which wrongly omitted N220F because
one unusable record affected its whole-mutation quality label.

There are **11 qualifying mutations**. In addition to the four above, these
are N270W, W325Q, C38S, E406Q, H223A, W325M, and E340A. The broader comparison
retains all **33 mutations with four usable records**, using their original
summaries, and gives more mixed results. Every selection decision is exported.

[View all 11 selected mutations](docs/figures/all_selected.png).

Records and contributor groups are not confirmed independent biological
replicates. N220F records 616 and 693 share the same kcat and KM values but
have different WT-reference information. Whether they share an experiment is
unresolved. The ten-record own-WT-reference check excludes 616 and retains
the same direction of change.

Excluded records are not described as proven noise. Missing information,
expression failure, and conflicting measurements are different issues.

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
highlights the most repeated qualifying mutations and `all_selected.png`
shows the entire selected set.

## Use the comparison code with other inputs

An invented example shows the required table format:

```bash
bglb-compare --consensus examples/synthetic/consensus.csv --measurements examples/synthetic/measurements.csv --predictions examples/synthetic/predictions.csv --policy selection_policy.json --model demo_model --dataset-label SYNTHETIC-DEMO --out demo_results
```

This example provides no evidence about CatPred. The generic command accepts
other audited tables in the same format. It checks record counts, grouping,
and measured summaries before comparing them with model results. The
[input guide](docs/methods.md#required-inputs) describes the fields.

## Scope and availability

See [the methods](docs/methods.md) for exact rules and
[the data guide](data/bglb/README.md) for sources and field definitions.
This release reproduces the comparison from previously audited inputs.
It does not run CatPred, train models, or recreate the original raw-data audit.

Dataset dates, CatPred settings, and the public-source check are documented in
the methods. The CatPred results come from a May 2026 export; its exact code
revision is unknown. This analysis is therefore not a test of CatPred's
current release or of the complete current D2D database.

The model results were generated with
[CatPred](https://github.com/maranasgroup/CatPred), developed by its authors.
This is not an official D2D or CatPred repository.

D2D provides the public experimental measurements. This project checks and
compares those measurements with CatPred results. On September 24, 2026, all
1,519 records in our dataset were found through D2D's public data endpoint,
and all kinetic values matched. All 86 entries for the 11 selected mutations,
including the 63 retained measurements, still have D2D approval and match the
public records in every compared field. Some other records have changed
approval flags or temperature measurements; the
[source-check results](docs/d2d_source_check.json) record those differences.

The release includes numerical measurements and audit decisions for all
1,519 source records, 968 usable mutant measurements, summaries for 756
mutations, and CatPred results for the 615 mutations with usable targets.
Contributor names are replaced with group IDs. Names, institutions, lab
notebooks, and instrument files are not included. Each source record keeps its
D2D ID so the numerical inputs can be checked against the public database.

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
