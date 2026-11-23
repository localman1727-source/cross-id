# -*- coding: utf-8 -*-
"""
CROSS-ID — one command from raw LC-MS/MS data to one trusted label per feature.

    python crossid.py init   "<raw data folder>" [--name <study>]   # write config/<study>.json
    python crossid.py status --config config/<study>.json          # which steps are done
    python crossid.py run    --config config/<study>.json          # run every step that is not done yet

Steps (in order)
  1 stage1       mzmine peak picking + MS-DIAL + candidate lists + merge      (run_stage1.py)
  2 run_dates    check whether acquisition date splits the samples            (check_run_date_effect.py)
  3 compare      sample comparison and selection categories                   (compare_samples.py)
  4 targets      pick annotation targets → targets.mgf (shared engine input)  (select_annotation_targets.py)
  5 sirius       SIRIUS formulas → fingerprints → structures (BIO only)       EXTERNAL
  6 ms2query     MS2Query on consensus spectra                                (run_ms2query.py)
  7 msfinder_in  MS-FINDER input files                                        (prepare_msfinder_input.py)
  8 msfinder     MS-FINDER structure prediction                               EXTERNAL
  9 msfinder_sum MS-FINDER summary                                            (summarize_msfinder.py)
 10 grade        engine agreement grade A–E                                   (grade_target_annotations.py --v2)
 11 occurrence   LOTUS occurrence + natural-product source flags              (flag_np_source.py)
 12 labels       one trusted label table per study                            (this file)

Why one runner?
  Each step reads files written by the step before. Run them out of order, or skip one, and the next step
  quietly reads an old file and still finishes without an error — the result is wrong and nothing says so.
  The runner checks each step's output before moving on and resumes from the first unfinished step.

EXTERNAL steps (SIRIUS needs an account login; MS-FINDER is a separate program)
  If the program is found and --run-external is given, the runner starts it and waits.
  Otherwise it stops, prints the exact command, and `run` picks up from there next time.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "cross_id"
sys.path.insert(0, str(SCRIPTS))
from study_config import load_study, utf8_console  # noqa: E402

# External tool locations (Windows defaults; override by environment variable or option)
SIRIUS_EXE = Path(os.environ.get("CROSSID_SIRIUS", r"C:\SIRIUS\sirius\sirius.exe"))
SIRIUS_WORK = Path(os.environ.get("CROSSID_SIRIUS_WORK", r"C:\SIRIUS\work"))
MSFINDER_EXE = Path(os.environ.get("CROSSID_MSFINDER", r"C:\MSFINDER\console_build\MsfinderConsoleApp.exe"))
MSFINDER_WORK = Path(os.environ.get("CROSSID_MSFINDER_WORK", r"C:\MSFINDER\work"))
MS2QUERY_PY = os.environ.get("CROSSID_MS2QUERY_PYTHON",
                             str(Path(os.environ.get("LOCALAPPDATA", "")) / "venvs" / "ms2query" / "Scripts" / "python.exe"))
TAG = "targets_consensus"


# ── step definitions ──────────────────────────────────────────────────────────
def _any(run, pattern):
    return bool(glob.glob(str(run / pattern)))


def msfinder_work(cfg):
    return MSFINDER_WORK / cfg["study_id"] / TAG


def msfinder_done(cfg):
    run = cfg["run_dir"]
    return (_any(run, f"msfinder/{TAG}/raw_results/Structure result-*.txt")
            or bool(glob.glob(str(msfinder_work(cfg) / "out_full" / "Structure result-*.txt"))))


def steps(cfg, conf):
    py = sys.executable
    run = cfg["run_dir"]
    S = lambda name: str(SCRIPTS / name)  # noqa: E731
    return [
        dict(key="stage1", title="mzmine peak picking + MS-DIAL + merge",
             cmd=[py, S("run_stage1.py"), *conf, "--run"],
             done=lambda: _any(run, "results/*_merged.csv")),
        dict(key="run_dates", title="acquisition-date check",
             cmd=[py, S("check_run_date_effect.py"), *conf],
             done=lambda: _any(run, "compare/*_run_date_check.json")),
        dict(key="compare", title="sample comparison and selection categories",
             cmd=[py, S("compare_samples.py"), *conf],
             done=lambda: _any(run, "compare/*_selection_categories.csv")),
        dict(key="targets", title="annotation targets (shared engine input)",
             cmd=[py, S("select_annotation_targets.py"), *conf],
             done=lambda: _any(run, "targets/*_targets.mgf")),
        dict(key="sirius", title="SIRIUS (external)", external=True,
             cmd=sirius_cmd(cfg), runner=lambda: run_sirius(cfg),
             done=lambda: _any(run, "sirius/summary/structure_identifications.tsv")),
        dict(key="ms2query", title="MS2Query (consensus spectra)",
             cmd=[MS2QUERY_PY, S("run_ms2query.py"), *conf, "--targets", "--spectrum", "consensus"],
             done=lambda: _any(run, f"ms2query/*_ms2query_{TAG}_results.csv")),
        dict(key="msfinder_in", title="MS-FINDER input",
             cmd=[py, S("prepare_msfinder_input.py"), *conf, "--targets", "--spectrum", "consensus"],
             done=lambda: (msfinder_work(cfg) / "msfinder_method.ini").exists() or msfinder_done(cfg)),
        dict(key="msfinder", title="MS-FINDER (external)", external=True,
             cmd=msfinder_cmd(cfg), runner=lambda: run_msfinder(cfg),
             done=lambda: msfinder_done(cfg)),
        dict(key="msfinder_sum", title="MS-FINDER summary",
             cmd=[py, S("summarize_msfinder.py"), *conf, "--tag", TAG],
             done=lambda: _any(run, f"msfinder/{TAG}/*_msfinder_summary.csv")),
        dict(key="grade", title="engine agreement grade A–E",
             cmd=[py, S("grade_target_annotations.py"), *conf, "--v2"],
             done=lambda: _any(run, "targets/*_targets_graded_v2.csv")),
        dict(key="occurrence", title="LOTUS occurrence and source flags",
             cmd=[py, S("flag_np_source.py"), *conf],
             done=lambda: _any(run, "targets/*_targets_graded_v2_npsource.csv")),
        dict(key="labels", title="trusted label table", runner=lambda: write_labels(cfg),
             cmd=["(built in)"],
             done=lambda: _any(run, "labels/*_cross_id_labels.csv")),
    ]


# ── external tools ────────────────────────────────────────────────────────────
def sirius_cmd(cfg):
    w = SIRIUS_WORK / cfg["study_id"]
    return [str(SIRIUS_EXE), "--cores=12", "-i", str(w / "targets.mgf"), "-o", str(w / "targets.sirius"),
            "formulas", "--profile", "qtof", "fingerprints", "classes", "structures", "--db=BIO", "--exp=OFF",
            "summaries", "--output", str(w / "summary")]


def run_sirius(cfg):
    """BIO structure database only, PubChem expansion off — so SIRIUS searches the same kind of
    structures as MS-FINDER and engine agreement stays fair."""
    w = SIRIUS_WORK / cfg["study_id"]
    w.mkdir(parents=True, exist_ok=True)
    mgf = sorted(glob.glob(str(cfg["run_dir"] / "targets" / "*_targets.mgf")))[0]
    shutil.copy(mgf, w / "targets.mgf")
    env = dict(os.environ, JAVA_TOOL_OPTIONS="-Xmx8g")
    log = w / "run.log"
    with log.open("w", encoding="utf-8") as fh:
        subprocess.run(sirius_cmd(cfg), stdout=fh, stderr=subprocess.STDOUT, env=env)
    text = log.read_text(encoding="utf-8", errors="replace")
    # exit code 0 is not proof of success: SIRIUS prints these and still exits 0
    bad = [b for b in ("Unmatched argument", "SEVERE", "403") if b in text]
    if bad or "Workflow has been finished" not in text:
        print(f"   SIRIUS log looks wrong ({', '.join(bad) or 'no finish line'}) — see {log}")
        return False
    dst = cfg["run_dir"] / "sirius" / "summary"
    dst.mkdir(parents=True, exist_ok=True)
    for f in (w / "summary").glob("*.tsv"):
        shutil.copy(f, dst / f.name)
    return True


def collect_external(cfg, key):
    """If the user ran SIRIUS / MS-FINDER by hand, copy their results next to the study."""
    if key == "sirius":
        src = SIRIUS_WORK / cfg["study_id"] / "summary"
        if (src / "structure_identifications.tsv").exists():
            dst = cfg["run_dir"] / "sirius" / "summary"
            dst.mkdir(parents=True, exist_ok=True)
            for f in src.glob("*.tsv"):
                shutil.copy(f, dst / f.name)
            print(f"   collected SIRIUS summaries from {src}")
    if key == "msfinder":
        src = msfinder_work(cfg) / "out_full"
        files = list(src.glob("*result-*.txt"))
        if files:
            dst = cfg["run_dir"] / "msfinder" / TAG / "raw_results"
            dst.mkdir(parents=True, exist_ok=True)
            for f in files:
                shutil.copy(f, dst / f.name)


def msfinder_cmd(cfg):
    w = msfinder_work(cfg)
    return [str(MSFINDER_EXE), "predict", "-i", str(w / "input"), "-o", str(w / "out_full"),
            "-m", str(w / "msfinder_method.ini")]


def run_msfinder(cfg):
    subprocess.run(msfinder_cmd(cfg))
    # keep a copy of the raw results next to the study so another PC can continue
    dst = cfg["run_dir"] / "msfinder" / TAG / "raw_results"
    dst.mkdir(parents=True, exist_ok=True)
    for f in (msfinder_work(cfg) / "out_full").glob("*result-*.txt"):
        shutil.copy(f, dst / f.name)
    return msfinder_done(cfg)


# ── final table ───────────────────────────────────────────────────────────────
def write_labels(cfg):
    """One trusted label per feature: the graded table with occurrence flags, under a stable name."""
    run = cfg["run_dir"]
    src = sorted(glob.glob(str(run / "targets" / "*_targets_graded_v2_npsource.csv")))[0]
    out = run / "labels"
    out.mkdir(exist_ok=True)
    dst = out / f"{cfg['study_id']}_cross_id_labels.csv"
    shutil.copy(src, dst)
    (out / "README.txt").write_text(
        "One row per target feature.\n"
        "grade (A–E) = how many independent engines agree on the structure;\n"
        "occurrence columns = whether that structure is reported in the same species, genus or family (LOTUS).\n"
        "Grade and occurrence are kept side by side, never merged into one score.\n"
        f"Source: {Path(src).name}\n", encoding="utf-8")
    print(f"   wrote {dst}")
    return True


# ── commands ──────────────────────────────────────────────────────────────────
def show(cmd):
    return " ".join(f'"{c}"' if " " in str(c) else str(c) for c in cmd)


def cmd_status(cfg, conf):
    print(f"\nCROSS-ID · {cfg['study_id']} ({cfg['polarity']})  ·  {cfg['run_dir']}\n")
    first = None
    for i, s in enumerate(steps(cfg, conf), 1):
        ok = s["done"]()
        if not ok and first is None:
            first = s["key"]
        print(f"  {'✔' if ok else '·'} {i:2d} {s['key']:<13s} {s['title']}{'   ← EXTERNAL' if s.get('external') else ''}")
    print(f"\n  next: {first}" if first else "\n  all steps done")


def cmd_run(cfg, conf, args):
    keys = [s["key"] for s in steps(cfg, conf)]
    start = keys.index(args.start) if args.start else 0
    stop = keys.index(args.stop) if args.stop else len(keys) - 1
    for s in steps(cfg, conf)[start:stop + 1]:
        if args.only and s["key"] != args.only:
            continue
        if s.get("external") and not args.dry_run:
            collect_external(cfg, s["key"])
        if s["done"]() and not args.force:
            print(f"✔ {s['key']:<13s} already done — skipped")
            continue
        print(f"▶ {s['key']:<13s} {s['title']}\n   {show(s['cmd'])}")
        if args.dry_run:
            continue
        if s.get("external"):
            exe = Path(s["cmd"][0])
            if not (args.run_external and exe.exists()):
                print(f"\n■ Stopped at the external step '{s['key']}'.\n"
                      f"  Run the command above (or add --run-external if {exe.name} is installed),\n"
                      f"  then run `crossid.py run` again — it continues from here.")
                return 2
            if not s["runner"]():
                print(f"✘ {s['key']} did not finish — fix and rerun.")
                return 1
            continue
        ok = s["runner"]() if "runner" in s else (subprocess.run([str(c) for c in s["cmd"]]).returncode == 0)
        if not ok or not s["done"]():
            print(f"✘ {s['key']} did not produce its output — stopped here so later steps do not read stale files.")
            return 1
    print("\nDone." if not args.dry_run else "\n(dry run — nothing executed)")
    return 0


def main():
    utf8_console()
    ap = argparse.ArgumentParser(description="CROSS-ID: consensus-ranked, occurrence-supported metabolite labels")
    sub = ap.add_subparsers(dest="command", required=True)
    p_init = sub.add_parser("init", help="create a study config from a raw-data folder")
    p_init.add_argument("raw_dir")
    p_init.add_argument("--name", default=None)
    for name in ("status", "run"):
        p = sub.add_parser(name)
        p.add_argument("--config", required=True)
        if name == "run":
            p.add_argument("--from", dest="start", default=None, help="start at this step")
            p.add_argument("--to", dest="stop", default=None, help="stop after this step")
            p.add_argument("--only", default=None, help="run only this step")
            p.add_argument("--force", action="store_true", help="rerun steps even if their output exists")
            p.add_argument("--run-external", action="store_true", help="start SIRIUS / MS-FINDER if installed")
            p.add_argument("--dry-run", action="store_true", help="print the commands only")
    args = ap.parse_args()

    if args.command == "init":
        extra = ["--name", args.name] if args.name else []
        return subprocess.run([sys.executable, str(SCRIPTS / "init_study.py"), args.raw_dir, *extra]).returncode

    cfg = load_study(args.config)
    conf = ["--config", str(Path(args.config).resolve())]
    if args.command == "status":
        cmd_status(cfg, conf)
        return 0
    return cmd_run(cfg, conf, args)


if __name__ == "__main__":
    sys.exit(main())
