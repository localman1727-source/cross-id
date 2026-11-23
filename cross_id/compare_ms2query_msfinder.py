# -*- coding: utf-8 -*-
"""
MS2Query(한 런 무늬 vs 합본 무늬) 비교 + MS-FINDER 와 같은 뼈대를 고른 봉우리 추리기
=====================================================================================
두 질문에 답한다.
  ① 합본 무늬(ms2_consensus.py)가 한 런 무늬보다 이름을 더 잘 붙이나?
     → 봉우리마다 두 방식의 1순위 점수(ms2query_model_prediction)를 나란히 놓고 오르내림을 센다.
       '같은 물질 후보'(전구 이온 질량 차 ≤ SAME_MZ_TOL Da)가 몇 개 늘었는지도 센다.
  ② 원리가 다른 두 도구가 같은 답을 냈나?
     → MS2Query 1순위의 뼈대(InChIKey 앞 14자)가 MS-FINDER 1순위와 같은지, 상위 10 안에 있는지.
       MS2Query 는 '실측 스펙트럼 사전에서 닮은 것 찾기', MS-FINDER 는 '구조 목록을 조각으로 쪼개 맞추기'라
       둘이 같은 뼈대를 고르면 서로 독립적인 근거 두 개가 된다.

'왜' 점수 자르는 선을 정하지 않나?
  MS2Query 논문은 0.7 근처를 믿을 만한 친척의 출발선으로 쓰지만 이 자료로 검증한 값이 아니다.
  그래서 표에는 점수를 그대로 두고, 요약(json)에는 몇 개 구간별 개수만 남긴다.
⚠️ 이름은 전부 후보다. 질량 차가 크면 '그만큼 무엇이 붙거나 떨어진 친척'이지 같은 물질이 아니다.

입력  runs/<연구>/ms2query/<연구>_ms2query_{last,consensus}_results.csv · _consensus_input_index.csv
      runs/<연구>/msfinder/<연구>_msfinder_summary.csv · C:/MSFINDER/work/<연구>/out_full/Structure result-*.txt
산출  runs/<연구>/ms2query/<연구>_ms2query_comparison.csv · .json · .md
실행  py scripts/compare_ms2query_msfinder.py --config config/<연구>.json
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
MSF_WORK = Path(r"C:\MSFINDER\work")
SAME_MZ_TOL = 0.01      # Da. 이 안이면 '같은 물질 후보'(분자량이 같다)
BANDS = [0.9, 0.8, 0.7, 0.6]


def top1(csv):
    d = pd.read_csv(csv, dtype={"feature_id": str})
    d = d[d["inchikey"].notna() | d["ms2query_model_prediction"].notna()]
    return d.groupby("feature_id", sort=False).head(1).set_index("feature_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_study(args.config)
    run, sid = cfg["run_dir"], cfg["study_id"]
    qd = run / "ms2query"

    L, C = top1(qd / f"{sid}_ms2query_last_results.csv"), top1(qd / f"{sid}_ms2query_consensus_results.csv")
    keep = ["ms2query_model_prediction", "precursor_mz_difference", "inchikey", "analog_compound_name", "smiles",
            "cf_superclass", "cf_class", "npc_pathway_results", "npc_superclass_results", "npc_class_results"]
    idx = pd.read_csv(qd / f"{sid}_ms2query_consensus_input_index.csv", dtype={"feature_id": str}).set_index("feature_id")
    out = idx[["precursor_mz", "rt_min", "mzmine_adduct", "n_spectra", "n_used", "median_cosine", "mode"]].rename(
        columns={"mode": "consensus_mode", "n_used": "consensus_n_used"})
    out = out.join(L[keep].add_prefix("last_")).join(C[keep].add_prefix("cons_"))
    out["score_change_cons_minus_last"] = out["cons_ms2query_model_prediction"] - out["last_ms2query_model_prediction"]
    sk = lambda s: s.fillna("").str[:14]
    out["last_cons_same_skeleton"] = (sk(out["last_inchikey"]) == sk(out["cons_inchikey"])) & sk(out["cons_inchikey"]).ne("")
    for m in ("last", "cons"):
        out[f"{m}_same_mass"] = out[f"{m}_precursor_mz_difference"].abs() <= SAME_MZ_TOL

    # MS-FINDER: 1순위 + 상위 10 뼈대
    msf = pd.read_csv(run / "msfinder" / f"{sid}_msfinder_summary.csv", dtype={"feature_id": str}).set_index("feature_id")
    sfile = sorted((MSF_WORK / sid / "out_full").glob("Structure result-*.txt"), key=lambda p: p.stat().st_mtime)[-1]
    sr = pd.read_csv(sfile, sep="\t", dtype=str, keep_default_na=False, quoting=3)
    top10 = {fid: [k[:14] for k in g["InChIKey"].head(10) if k] for fid, g in sr.groupby("File name", sort=False)}
    out = out.join(msf[["msf_top_name", "msf_top_inchikey", "msf_top_formula", "msf_top_in_relatives",
                        "mzmine_name", "mzmine_inchikey", "msdial_suggestion"]])
    for m in ("last", "cons"):
        k = sk(out[f"{m}_inchikey"])
        out[f"{m}_agrees_msf_top1"] = k.ne("") & (k == sk(out["msf_top_inchikey"]))
        out[f"{m}_in_msf_top10"] = [bool(a) and a in top10.get(f, []) for f, a in zip(out.index, k)]
    out = out.reset_index().rename(columns={"index": "feature_id"})
    out.to_csv(qd / f"{sid}_ms2query_comparison.csv", index=False, encoding="utf-8-sig")

    n = len(out)
    def bands(col):
        s = out[col]
        return {f">={b}": int((s >= b).sum()) for b in BANDS}
    cnt = {
        "features": n,
        "consensus_mode_counts": out["consensus_mode"].value_counts().to_dict(),
        "score_bands_last": bands("last_ms2query_model_prediction"),
        "score_bands_consensus": bands("cons_ms2query_model_prediction"),
        "same_mass_and_score>=0.7_last": int((out.last_same_mass & (out.last_ms2query_model_prediction >= 0.7)).sum()),
        "same_mass_and_score>=0.7_consensus": int((out.cons_same_mass & (out.cons_ms2query_model_prediction >= 0.7)).sum()),
        "consensus_score_higher": int((out.score_change_cons_minus_last > 0).sum()),
        "consensus_score_lower": int((out.score_change_cons_minus_last < 0).sum()),
        "consensus_score_equal": int((out.score_change_cons_minus_last == 0).sum()),
        "median_score_change": float(out.score_change_cons_minus_last.median()),
        "last_and_consensus_same_skeleton": int(out.last_cons_same_skeleton.sum()),
        "agrees_msfinder_top1_last": int(out.last_agrees_msf_top1.sum()),
        "agrees_msfinder_top1_consensus": int(out.cons_agrees_msf_top1.sum()),
        "in_msfinder_top10_last": int(out.last_in_msf_top10.sum()),
        "in_msfinder_top10_consensus": int(out.cons_in_msf_top10.sum()),
    }
    # 합본이 실제로 쓰인 봉우리(무늬가 엇갈리던 곳)에서만 따로
    sub = out[out.consensus_mode == "consensus"]
    cnt["among_consensus_built"] = {
        "n": int(len(sub)),
        "consensus_score_higher": int((sub.score_change_cons_minus_last > 0).sum()),
        "consensus_score_lower": int((sub.score_change_cons_minus_last < 0).sum()),
        "median_score_change": float(sub.score_change_cons_minus_last.median())}
    info = {"same_mz_tol_Da": SAME_MZ_TOL, "msfinder_structure_file": str(sfile), "counts": cnt}
    (qd / f"{sid}_ms2query_comparison.json").write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")

    agree = out[out.cons_agrees_msf_top1 | out.last_agrees_msf_top1].sort_values("cons_ms2query_model_prediction", ascending=False)
    md = [f"# {sid} — MS2Query(한 런 vs 합본) · MS-FINDER 대조", "",
          "이름은 전부 후보다. 질량 차가 0 근처면 같은 물질 후보, 크면 친척.", "",
          "| 항목 | 값 |", "|---|---|"] + [f"| {k} | {v} |" for k, v in cnt.items()] + ["",
          f"## MS2Query 와 MS-FINDER 1순위 뼈대가 같은 봉우리 ({len(agree)})", "",
          "| feature | m/z | RT | 이름(MS2Query 합본) | 점수 한런→합본 | 질량차 | MS-FINDER 1순위 | 근연 목록 |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in agree.iterrows():
        md.append(f"| {r.feature_id} | {r.precursor_mz} | {r.rt_min} | {r.cons_analog_compound_name} | "
                  f"{r.last_ms2query_model_prediction}→{r.cons_ms2query_model_prediction} | {r.cons_precursor_mz_difference} | "
                  f"{r.msf_top_name} | {r.msf_top_in_relatives} |")
    (qd / f"{sid}_ms2query_comparison.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(cnt, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
