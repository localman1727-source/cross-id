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

## Requirements

- Python 3.13 (tested on Windows; pandas, numpy, rdkit — see each script header)
- External tools installed separately, each under its own licence: mzmine 3/4, MS-DIAL, MS2Query, MS-FINDER, SIRIUS (account login required).
- Developed and tested on Windows with Waters SYNAPT XS (MSe) data. Default tool paths are Windows paths and can be changed in the scripts or via arguments.

No raw data, study data or engine binaries are included in this repository.

## Citation

If you use CROSS-ID, please cite this repository. A methods paper is in preparation.

## Licence

MIT © 2026 Jimin Hyun
