# -*- coding: utf-8 -*-
"""
MS-FINDER 결과 정리 — 봉우리마다 1순위 분자식·구조 후보와, 그 후보가 근연 식물에서 보고된 것인지
=======================================================================================================
'왜' 1순위만 믿으면 안 되나?
  MS-FINDER 는 조각무늬가 후보 구조에서 얼마나 설명되는지로 순위를 매긴다. 같은 분자식의 이성질체
  (예: kaempferol 과 luteolin)는 조각이 거의 같아 점수 차가 작다. 그래서 1순위 이름과 함께
  ① 2순위와의 점수 차 ② 상위 10개 중 1순위와 분자식이 같은 후보 수를 남겨, 이름을 읽는 사람이
  '같은 분자식의 다른 구조일 수 있음'을 볼 수 있게 한다. 점수 차를 어디서 자를지는 정하지 않는다.

'왜' 근연 식물 목록과 대조하나?
  MS-FINDER 의 구조 목록은 모든 생물의 천연물을 담고 있어 해양·미생물 물질도 후보에 오른다.
  LOTUS 에서 만든 근연 속·과 후보 목록(build_candidate_db.py)에 같은 뼈대(InChIKey 앞 14자)가
  있으면 '가까운 식물에서 실제 보고된 구조'라는 보조 근거가 된다. 없다고 틀린 것은 아니다.

⚠️ 이 표의 이름은 전부 후보(조각무늬 + 구조 목록 기반 예측)이지 표준품으로 확인한 동정이 아니다.
⚠️ 첨가 이온을 [M+H]+ 로 고정해 돌렸다(prepare_msfinder_input.py). mzmine 이 다른 첨가 이온으로
   본 봉우리는 mzmine_adduct 칸에 남아 있다 — 그런 행의 분자식은 따로 확인이 필요하다.

입력  C:/MSFINDER/work/<연구>/out_full/{Formula,Structure} result-*.txt  (가장 최근 것)
      runs/<연구>/msfinder/input_index.csv · runs/<연구>/sirius/<연구>_sirius_selection.csv
      runs/<연구>/results/<연구>_merged.csv · runs/<연구>/db/<연구>_candidates_{relatives,genus,family}.csv
산출  runs/<연구>/msfinder/<연구>_msfinder_summary.csv · _summary.json · _summary.md
실행  py scripts/summarize_msfinder.py --config config/<연구>.json [--out-subdir out_full] [--tag targets_consensus]
      --tag 를 주면 작업 폴더 C:/MSFINDER/work/<연구>/<tag>/, 기록 폴더 runs/<연구>/msfinder/<tag>/ 를 쓰고
      신호 크기·QC 칸은 선발 목록(runs/<연구>/targets/<연구>_targets.csv)에서 가져온다(범주 칸 포함).
"""
import sys
import json
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console

utf8_console()
WORK_ROOT = Path(r"C:\MSFINDER\work")


def latest(folder, prefix):
    fs = sorted(Path(folder).glob(f"{prefix}*.txt"), key=lambda p: p.stat().st_mtime)
    if not fs:
        raise SystemExit(f"[중단] {folder} 에 {prefix}*.txt 없음")
    return fs[-1]


def read_tsv(p):
    # 구조 결과의 머리줄은 'Structure rank' 자리에 실제로 이름이 들어 있다(v3.73 콘솔 출력 그대로).
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, quoting=3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--out-subdir", default="out_full")
    ap.add_argument("--tag", default=None, help="prepare_msfinder_input.py 의 tag (예: targets_consensus)")
    args = ap.parse_args()

    cfg = load_study(args.config)
    run, sid = cfg["run_dir"], cfg["study_id"]
    fixed = "[M+H]+" if cfg["polarity"].upper() == "POS" else "[M-H]-"   # prepare_msfinder_input.py 와 같은 규칙
    work_out = (WORK_ROOT / sid / args.tag if args.tag else WORK_ROOT / sid) / args.out_subdir
    mdir = run / "msfinder" / args.tag if args.tag else run / "msfinder"
    fpath, spath = latest(work_out, "Formula result-"), latest(work_out, "Structure result-")
    fr, sr = read_tsv(fpath), read_tsv(spath)
    sr = sr.rename(columns={"Structure rank": "name"})

    # 근연 식물 목록: 뼈대(앞 14자) → 어느 목록에 있었나
    rel = {}
    for tag in ("relatives", "genus", "family"):
        p = run / "db" / f"{sid}_candidates_{tag}.csv"
        if p.exists():
            d = pd.read_csv(p, dtype=str, keep_default_na=False)
            for k, n, sp in zip(d["inchikey"], d["name"], d.get("source_species", [""] * len(d))):
                if k:
                    rel.setdefault(k[:14], (tag, n, sp))

    idx = pd.read_csv(mdir / "input_index.csv", dtype={"feature_id": str})
    if args.tag and args.tag.startswith("targets"):
        sel = pd.read_csv(run / "targets" / f"{sid}_targets.csv", dtype={"feature_id": str}).rename(
            columns={"feature_id": "id"})[["id", "height", "qc_rsd", "category"]]
    else:
        sel = pd.read_csv(run / "sirius" / f"{sid}_sirius_selection.csv", dtype={"id": str})
    mg = pd.read_csv(run / "results" / f"{sid}_merged.csv", dtype={"mzmine_id": str}, low_memory=False)
    mg = mg[mg["mzmine_id"].notna()].drop_duplicates("mzmine_id")[
        ["mzmine_id", "mzmine_name", "mzmine_formula", "mzmine_inchikey", "msdial_suggestion"]]

    rows = []
    for fid, g in sr.groupby("File name", sort=False):
        g = g.reset_index(drop=True)
        g["score"] = pd.to_numeric(g["Total score"], errors="coerce")
        top = g.iloc[0]
        top10 = g.head(10)
        rel_hits = [(i + 1, rel[k[:14]]) for i, k in enumerate(top10["InChIKey"]) if k and k[:14] in rel]
        rows.append({
            "feature_id": fid,
            "msf_top_name": top["name"], "msf_top_formula": top["Formula"],
            "msf_top_score": top["score"], "msf_top_inchikey": top["InChIKey"],
            "msf_top_ontology": top["Ontology"], "msf_top_smiles": top["SMILES"],
            "msf_score_gap_to_2nd": (top["score"] - g.iloc[1]["score"]) if len(g) > 1 else None,
            "msf_n_same_formula_top10": int((top10["Formula"] == top["Formula"]).sum()),
            "msf_n_candidates": len(g),
            "msf_top_in_relatives": bool(top["InChIKey"] and top["InChIKey"][:14] in rel),
            "msf_best_relatives_rank": rel_hits[0][0] if rel_hits else None,
            "msf_best_relatives_name": rel_hits[0][1][1] if rel_hits else None,
            "msf_best_relatives_list": rel_hits[0][1][0] if rel_hits else None,
            "msf_best_relatives_species": rel_hits[0][1][2] if rel_hits else None,
            "msf_top_databases": top["Databases"],
        })
    st = pd.DataFrame(rows)

    # 분자식 1순위(구조 후보가 없어도 분자식은 있을 수 있다)
    fr1 = fr.groupby("File name", sort=False).head(1)[["File name", "Formula", "Formula score", "Mass error"]]
    fr1.columns = ["feature_id", "msf_formula_top", "msf_formula_score", "msf_formula_mass_error_Da"]

    out = (idx.merge(sel.rename(columns={"id": "feature_id"})[[c for c in ("feature_id", "height", "qc_rsd", "category") if c in sel.columns or c == "feature_id"]], on="feature_id", how="left")
              .merge(fr1, on="feature_id", how="left")
              .merge(st, on="feature_id", how="left")
              .merge(mg.rename(columns={"mzmine_id": "feature_id"}), on="feature_id", how="left"))
    out["adduct_differs_from_mzmine"] = out["mzmine_adduct"].fillna("").ne("") & out["mzmine_adduct"].ne(fixed)
    out = out.sort_values("height", ascending=False)

    out.to_csv(mdir / f"{sid}_msfinder_summary.csv", index=False, encoding="utf-8-sig")

    n = len(out)
    counts = {
        "features_input": int(n),
        "with_formula": int(out["msf_formula_top"].notna().sum()),
        "with_structure_candidate": int(out["msf_top_name"].notna().sum()),
        "top1_in_relatives_list": int(out["msf_top_in_relatives"].fillna(False).astype(bool).sum()),
        "any_top10_in_relatives_list": int(out["msf_best_relatives_rank"].notna().sum()),
        "top1_has_same_formula_rivals_in_top10": int((out["msf_n_same_formula_top10"] > 1).sum()),
        "adduct_differs_from_mzmine": int(out["adduct_differs_from_mzmine"].sum()),
        "mzmine_also_named": int(out["mzmine_name"].notna().sum()),
        "mzmine_named_and_same_skeleton_as_msf_top1": int(sum(
            1 for a, b in zip(out["mzmine_inchikey"].fillna(""), out["msf_top_inchikey"].fillna(""))
            if a and b and a[:14] == b[:14])),
    }
    info = {"formula_result": str(fpath), "structure_result": str(spath),
            "msfinder": "MsfinderConsoleApp (MsdialWorkbench master bb90e0e5, net48 자체 빌드) + MS-FINDER v3.73 Resources",
            "precursor_type": f"{fixed} fixed", "pubchem_online": False, "mines": False,
            "relatives_lists": sorted({v[0] for v in rel.values()}), "counts": counts}
    (mdir / f"{sid}_msfinder_summary.json").write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")

    md = [f"# {sid} — MS-FINDER 구조 후보 요약", "",
          f"모든 이름은 후보다(조각무늬 + 공개 구조 목록 기반 예측, 표준품 확인 아님). 첨가 이온은 {fixed} 로 고정했다.", "",
          "| 항목 | 개수 |", "|---|---|"] + [f"| {k} | {v} |" for k, v in counts.items()] + ["",
          "## 신호 큰 순 상위 30", "",
          "| feature | m/z | RT | 1순위 후보 | 분자식 | 점수 | 2순위와 차 | 같은 분자식(상위10) | 근연 목록 최고순위 | mzmine 이름 |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in out.head(30).iterrows():
        md.append(f"| {r.feature_id} | {r.precursor_mz} | {r.rt_min} | {r.msf_top_name} | {r.msf_top_formula} | "
                  f"{r.msf_top_score} | {r.msf_score_gap_to_2nd} | {r.msf_n_same_formula_top10} | "
                  f"{'' if pd.isna(r.msf_best_relatives_rank) else f'{int(r.msf_best_relatives_rank)}: {r.msf_best_relatives_name}'} | "
                  f"{'' if pd.isna(r.mzmine_name) else r.mzmine_name} |")
    (mdir / f"{sid}_msfinder_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(counts, ensure_ascii=False))
    print(f"[생성] {mdir}")


if __name__ == "__main__":
    main()
