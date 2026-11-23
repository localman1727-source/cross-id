# CROSS-ID

**Consensus-Ranked, Occurrence-Supported Scoring for IDentification**
A trust-based metabolite labelling workflow for untargeted LC-MS/MS of natural products and foods.

> Status: early private release (2026-11). Code comments are mostly in Korean; English documentation is in progress.

## What it does

Single annotation engines often disagree, and a top-ranked name can be chemically or biologically implausible.
CROSS-ID does not trust any single engine. It covers the whole path from raw data to a trusted label —
peak picking with mzmine plus four annotation engines — and for every feature it:

1. **Picks peaks with mzmine** — one peak list for all samples, from a batch generated out of the study config and a method template (`templates/`).
2. **Runs four annotation engines** on the same targets and the same consensus spectra — MS-DIAL, MS2Query, MS-FINDER and SIRIUS.
3. **Grades engine agreement A–E** — from "all engines agree on the structure" to "none agree".
4. **Adds occurrence evidence** — whether the structure has been reported in the same species, genus or family (LOTUS structure–organism records) and the source papers.
5. **Writes one trusted label per feature** — name, grade, evidence and the engines that agreed.

Only graded, evidenced labels move on to biological interpretation.

## Quick start

```
python crossid.py init   "D:\raw\my_study"                 # write config/<study>.json from the raw-data folder
python crossid.py status --config config/<study>.json       # which steps are done
python crossid.py run    --config config/<study>.json       # run every step that is not done yet
```

`run` goes through twelve steps in order — mzmine peak picking and MS-DIAL, acquisition-date check, sample
comparison, target selection, SIRIUS, MS2Query, MS-FINDER, engine-agreement grade, LOTUS occurrence, and the
final label table (`runs/<study>/labels/<study>_cross_id_labels.csv`). It checks each step's output before moving
on, skips steps that are already done, and stops instead of letting a later step read a stale file.

SIRIUS (account login) and MS-FINDER run as separate programs. If they are installed, add `--run-external`;
otherwise `run` stops at that step, prints the exact command, and continues from there next time.
Other options: `--from <step>`, `--to <step>`, `--only <step>`, `--force`, `--dry-run`.
Tool locations can be set with `CROSSID_SIRIUS`, `CROSSID_MSFINDER`, `CROSSID_MS2QUERY_PYTHON` and the `*_WORK` variables.

## Layout

```
crossid.py         single entry point (init / status / run)
cross_id/          step scripts called by crossid.py (study_config.py reads one JSON per study)
config/            study_template.json — copy and fill in, or use `crossid.py init`
templates/         mzmine batch templates (Waters MSe examples; add your own)
```

## Components — all free for academic use

CROSS-ID adds no paid software. Every component below can be used free of charge for academic and
non-commercial research. Two of them ask for a free account, and both offer paid licences for commercial use.

| Component | Role in CROSS-ID | Cost and licence | Get it |
|---|---|---|---|
| mzmine 4 | Peak picking, one peak list for all samples | Free COMMUNITY edition for academic use (free user account). Source code MIT. Commercial use: mzmine PRO | [mzmine.github.io](https://mzmine.github.io/) · [releases](https://github.com/mzmine/mzmine/releases) · [getting started](https://mzmine.github.io/mzmine_documentation/getting_started.html) |
| MS-DIAL | Spectral library matching | Free download (open source, RIKEN / UC Davis) | [MS-DIAL](https://systemsomicslab.github.io/compms/msdial/main.html) · [source](https://github.com/systemsomicslab/MsdialWorkbench) |
| MS-FINDER | In-silico fragmentation against structure databases | Free download (open source, same project as MS-DIAL). CROSS-ID calls the console app built from the MsdialWorkbench source | [MS-FINDER](https://systemsomicslab.github.io/compms/msfinder/main.html) · [source](https://github.com/systemsomicslab/MsdialWorkbench) |
| MS2Query | Machine-learning library search for analogues | Free, Apache-2.0. Install: `pip install ms2query` (separate Python environment recommended). Pretrained models on Zenodo | [GitHub](https://github.com/iomega/ms2query) · models: [positive](https://zenodo.org/records/13348638) · [negative](https://zenodo.org/records/13348683) |
| SIRIUS 6 (CSI:FingerID, CANOPUS) | Formula, structure fingerprint and compound class prediction | Free for academic / non-commercial use; free account, academic licence granted from the institutional e-mail domain. Commercial use: Bright Giant | [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) · [account and licence](https://boecker-lab.github.io/docs.sirius.github.io/account-and-license) |
| LOTUS | Structure–organism occurrence records (evidence, not an engine) | Open data (CC0) | [LOTUS](https://lotus.naturalproducts.net/) · frozen snapshot used: [Zenodo 19360665](https://zenodo.org/records/19360665) |
| Python 3.13 + pandas, numpy, rdkit | Glue, grading, tables | Free, open source | [python.org](https://www.python.org/downloads/) · `pip install pandas numpy rdkit` |

Check each tool's own licence before commercial or patent-related work.

## Requirements

- Developed and tested on Windows with Waters SYNAPT XS (MSe) data. Default tool locations are Windows paths;
  change them with `CROSSID_SIRIUS`, `CROSSID_MSFINDER`, `CROSSID_MS2QUERY_PYTHON` and the `*_WORK` variables.
- No raw data, study data or engine binaries are included in this repository.

## Citation

If you use CROSS-ID, please cite this repository. A methods paper is in preparation.

## Licence

MIT © 2026 Jimin Hyun
