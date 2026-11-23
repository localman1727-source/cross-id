# -*- coding: utf-8 -*-
"""
선발 대상(산지 특이·우세 + 공통)의 이름 근거를 한 줄로 모아 믿음 등급을 매긴다
==============================================================================
근거 셋을 봉우리마다 나란히 둔다.
  MS-DIAL   공개 스펙트럼 사전과 조각무늬 직접 대조(merged.csv 의 msdial_name)
  MS2Query  공개 실측 사전에서 같은 물질/친척 찾기(합본 무늬, run_ms2query.py --targets --spectrum consensus)
  MS-FINDER 구조 목록을 조각으로 쪼개 맞추기(합본 무늬, prepare_msfinder_input.py --targets --spectrum consensus)

등급 (높은 것부터, 먼저 맞는 것 하나)
  A_library        MS-DIAL 이 조각무늬로 이름을 붙였다
  B_two_tools      MS2Query 1순위가 같은 물질 후보(질량차 ≤ SAME_MZ_TOL, 점수 ≥ SCORE_SAME)이고
                   그 뼈대(InChIKey 앞 14자)가 MS-FINDER 상위 10 안에 있다 — 원리가 다른 두 도구가 같은 답
  C_one_tool       MS2Query 같은 물질 후보만, 또는 MS-FINDER 1순위가 근연 식물 목록에 있음만
  D_class_only     이름은 못 믿지만 MS2Query 친척의 물질 계열은 참고 가능(점수 ≥ SCORE_CLASS)
  E_unknown        위 어느 것도 아님
'왜' 이 문턱인가? 0.7·0.6 은 MS2Query 논문에서 친척 판단에 쓰는 언저리 값이고 이 자료로 검증한 값이
  아니다. 문턱은 등급을 나누는 데만 쓰고, 표에는 점수를 그대로 둔다.
⚠️ 등급은 '근거가 몇 겹인가'이지 확정이 아니다. A 도 표준품 확인 전에는 후보다(Level 2).
⚠️ 사전 편향 이름(도핑 스테로이드 등)은 등급과 무관하게 의심 칸(name_flag)에 표시한다.

v2 (--v2, 2026-11-06 결정 — 음이온 이전 연구 부터. 첫 연구 v1 산출은 그대로 둔다)
  ① SIRIUS(CSI:FingerID 1순위 뼈대)를 넷째 도구로 넣고, 같은 뼈대를 가리키는 도구 수(n_tools_agree)를 센다.
     B_multi_tool = MS-DIAL 이름이 없고, 원리가 다른 도구 둘 이상이 같은 뼈대.
     도구별 '지목' 조건: MS-DIAL 이름 · MS2Query 같은 물질 후보(v1 과 같은 문턱) · MS-FINDER 상위 10 ·
     SIRIUS 1순위. 1순위가 늘 나오는 MS-FINDER·SIRIUS 는 혼자서는 C 가 되지 않는다(근연 목록에 있을 때만, v1 규칙과 같게).
  ② 출현 근거(occurrence)를 등급과 나란히 둔다(합치거나 강등하지 않는다):
     species · genus · relatives · family · lotus_other · none  (LOTUS 2026-04-13, 후보 목록 db/*_candidates_*.csv)
     '왜' 강등하지 않나? 대상 수종는 종 단위 보고가 0 이라 강등하면 새 성분 후보까지 일괄로 밀린다.
     이름을 내걸 때만: 등급 C 이하 + occurrence none → name_policy = class_only(계열만 표시).
  ③ SIRIUS CANOPUS 계열(NPC pathway·class)은 칸으로만 남긴다 — D 등급 조건은 v1(MS2Query)과 같다.
  산출: runs/<연구>/targets/<연구>_targets_graded_v2.{csv,json,md}
  SIRIUS 요약 폴더는 --sirius-dir 로(기본 runs/<연구>/sirius/summary). 요약의 mappingFeatureId 가 대상
  FEATURE_ID 와 하나도 안 맞으면 멈춘다(조용히 SIRIUS 0건 되는 것 방지).

입력  runs/<연구>/targets/<연구>_targets.csv · runs/<연구>/results/<연구>_merged.csv
      runs/<연구>/ms2query/<연구>_ms2query_targets_consensus_results.csv
      runs/<연구>/msfinder/targets_consensus/<연구>_msfinder_summary.csv · C:/MSFINDER/work/<연구>/targets_consensus/out_full/Structure result-*.txt
산출  runs/<연구>/targets/<연구>_targets_graded.csv · .json · .md
실행  py scripts/grade_target_annotations.py --config config/<연구>.json
"""
import sys
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console

utf8_console()
MSF_WORK = Path(r"C:\MSFINDER\work")
SAME_MZ_TOL = 0.01
SCORE_SAME = 0.7
SCORE_CLASS = 0.6
# 공개 사전에 표준품이 몰려 높은 점수로 잘못 붙기 쉬운 이름(식물 시료에서 그대로 읽으면 안 됨)
SUSPECT_WORDS = ["nandrolone", "nortestosterone", "boldenone", "testosterone", "androst", "estradiol",
                 "dark web", "amoxicil", "ergotamin", "drug"]


OCC_ORDER = ["species", "genus", "relatives", "family", "lotus_other", "none"]


def occurrence_lookup(cfg, run, sid):
    """뼈대(InChIKey 앞 14자) → 가장 가까운 출현 표지. 후보 목록은 build_candidate_db.py 산출을 그대로 읽는다."""
    binom = " ".join(cfg["samples"][0]["scientific_name"].split()[:2])   # 대조는 두 낱말로만(PROCESS ③)
    occ = {}

    def put(keys, tag):
        for k in keys:
            if isinstance(k, str) and k and k not in occ:
                occ[k] = tag
    g = pd.read_csv(run / "db" / f"{sid}_candidates_genus.csv", dtype=str)
    put(g.loc[g.source_species.fillna("").str.contains(binom, regex=False), "inchikey"].str[:14], "species")
    for tag in ("genus", "relatives", "family"):
        f = run / "db" / f"{sid}_candidates_{tag}.csv"
        if f.exists():
            put(pd.read_csv(f, dtype=str, usecols=["inchikey"]).inchikey.str[:14], tag)
    lotus = Path.home() / ".npdb_cache" / "lotus_260413_frozen_metadata.csv.gz"
    if not lotus.exists():
        raise SystemExit(f"[중단] LOTUS 캐시 없음: {lotus} — build_candidate_db.py 를 먼저 돌리면 받아진다")
    lk = pd.read_csv(lotus, usecols=["structure_inchikey"], dtype=str).structure_inchikey.dropna().str[:14]
    put(lk.unique(), "lotus_other")
    return occ, binom


def load_sirius(sdir, ids):
    st = pd.read_csv(sdir / "structure_identifications.tsv", sep="\t", dtype=str, keep_default_na=False)
    if not st.mappingFeatureId.isin(ids).any():
        raise SystemExit("[중단] SIRIUS 요약의 mappingFeatureId 가 대상 FEATURE_ID 와 하나도 맞지 않습니다 — 입력 mgf 확인")
    st = st[st.structurePerIdRank == "1"].drop_duplicates("mappingFeatureId").set_index("mappingFeatureId")
    out = pd.DataFrame({"s_inchikey14": st.InChIkey2D.str[:14], "s_name": st.name, "s_formula": st.molecularFormula,
                        "s_adduct": st.adduct, "s_conf_exact": pd.to_numeric(st.ConfidenceScoreExact, errors="coerce"),
                        "s_conf_approx": pd.to_numeric(st.ConfidenceScoreApproximate, errors="coerce")})
    cf = sdir / "canopus_structure_summary.tsv"
    if cf.exists():
        c = pd.read_csv(cf, sep="\t", dtype=str, keep_default_na=False)
        c = c[c.formulaRank == "1"].drop_duplicates("mappingFeatureId").set_index("mappingFeatureId")
        out = out.join(pd.DataFrame({"s_npc_pathway": c["NPC#pathway"], "s_npc_class": c["NPC#class"]}), how="outer")
    return out


def grade_v2(cfg, run, sid, d, q_same, rel, top10, sfile, sirius_dir):
    sdir = Path(sirius_dir) if sirius_dir else run / "sirius" / "summary"
    if not (sdir / "structure_identifications.tsv").exists():
        raise SystemExit(f"[중단] SIRIUS 요약 없음: {sdir} — SIRIUS 를 먼저 돌리거나 --sirius-dir 지정")
    d = d.join(load_sirius(sdir, set(d.index)))
    occ, binom = occurrence_lookup(cfg, run, sid)
    # 2026-11-06 이전 연구: MS-FINDER 를 대상 일부에만 돌린 경우(msfinder_subset_features.csv) — 안 돌린 대상은
    # '지목 없음'이 아니라 '안 돌림'으로 칸을 남긴다(등급 계산은 같고, 읽을 때 구분).
    sub = run / "msfinder" / "targets_consensus" / "msfinder_subset_features.csv"
    if sub.exists():
        ran = set(pd.read_csv(sub, dtype={"feature_id": str}).feature_id)
        d["msfinder_ran"] = [f in ran for f in d.index]
    else:
        d["msfinder_ran"] = True

    def k14(x):
        return x[:14] if isinstance(x, str) and x else None

    rows = []
    for f, r in d.iterrows():
        picks = {}                                      # 도구 → 지목한 뼈대
        if isinstance(r.msdial_name, str) and k14(r.msdial_inchikey):
            picks["MS-DIAL"] = k14(r.msdial_inchikey)
        if bool(q_same.get(f, False)) and k14(r.q_inchikey):
            picks["MS2Query"] = k14(r.q_inchikey)
        if k14(r.s_inchikey14):
            picks["SIRIUS"] = k14(r.s_inchikey14)
        cands = set(picks.values()) | ({k14(r.msf_top_inchikey)} if k14(r.msf_top_inchikey) else set())
        best_k, best_tools = None, []
        for k in sorted(cands):
            tools = [t for t, v in picks.items() if v == k] + (["MS-FINDER"] if k in top10.get(f, set()) else [])
            if len(tools) > len(best_tools):
                best_k, best_tools = k, tools
        rows.append({"feature_id": f, "agree_inchikey14": best_k, "n_tools_agree": len(best_tools),
                     "tools_agree": "+".join(sorted(best_tools))})
    d = d.join(pd.DataFrame(rows).set_index("feature_id"))

    s_in_rel = d.s_inchikey14.map(lambda k: occ.get(k) in ("species", "genus", "relatives") if isinstance(k, str) else False)
    g = pd.Series("E_unknown", index=d.index)
    g[d.q_ms2query_model_prediction >= SCORE_CLASS] = "D_class_only"
    g[q_same | rel | s_in_rel] = "C_one_tool"
    g[d.n_tools_agree >= 2] = "B_multi_tool"
    g[d.msdial_name.notna()] = "A_library"
    d["grade_v2"] = g

    # 출현 근거는 등급을 매긴 뼈대 기준: 여러 도구 합의 > MS-DIAL > MS2Query > SIRIUS
    key = d.agree_inchikey14.where(d.n_tools_agree >= 2,
          d.msdial_inchikey.str[:14].where(d.msdial_name.notna(),
          d.q_inchikey.str[:14].where(q_same, d.s_inchikey14)))
    d["occ_inchikey14"] = key
    d["occurrence"] = key.map(lambda k: occ.get(k, "none") if isinstance(k, str) and k else "none")
    weak = d.grade_v2.isin(["C_one_tool", "D_class_only", "E_unknown"]) & (d.occurrence == "none")
    d["name_policy"] = np.where(weak, "class_only", "name")
    d.to_csv(run / "targets" / f"{sid}_targets_graded_v2.csv", encoding="utf-8-sig", index_label="feature_id")

    ct = pd.crosstab(d.grade_v2, d.occurrence).reindex(columns=[c for c in OCC_ORDER if c in set(d.occurrence)])
    info = {"version": "v2 (2026-11-06)", "species_binomial": binom, "sirius_dir": str(sdir),
            "thresholds": {"same_mz_tol_Da": SAME_MZ_TOL, "score_same": SCORE_SAME, "score_class": SCORE_CLASS},
            "grade_total": d.grade_v2.value_counts().to_dict(),
            "grade_by_category": {c: v.to_dict() for c, v in pd.crosstab(d.category, d.grade_v2).iterrows()},
            "grade_x_occurrence": {gname: row.to_dict() for gname, row in ct.iterrows()},
            "n_tools_agree": {int(k): int(v) for k, v in d.n_tools_agree.value_counts().sort_index().items()},
            "tools_agree_combos": d.loc[d.n_tools_agree >= 2, "tools_agree"].value_counts().to_dict(),
            "name_policy_class_only": int((d.name_policy == "class_only").sum()),
            "msfinder_ran": int(d.msfinder_ran.sum()), "msfinder_not_run": int((~d.msfinder_ran).sum()),
            "grade_by_msfinder_ran": {str(k): v.to_dict() for k, v in pd.crosstab(d.msfinder_ran, d.grade_v2).iterrows()},
            "msfinder_structure_file": str(sfile),
            "lotus": "frozen 2026-04-13, Zenodo 10.5281/zenodo.19360665"}
    (run / "targets" / f"{sid}_targets_graded_v2.json").write_text(
        json.dumps(info, ensure_ascii=False, indent=1, default=str), encoding="utf-8")

    md = [f"# {sid} — 선발 대상 이름 근거 등급 v2 (도구 일치 × 출현 근거)", "",
          "등급은 근거의 겹 수이지 확정이 아니다. 출현 근거는 등급과 따로 읽는다.", "",
          "## 등급 × 출현 근거", "", "| 등급 | " + " | ".join(ct.columns) + " |", "|---" * (len(ct.columns) + 1) + "|"] + \
         [f"| {i} | " + " | ".join(str(v) for v in ct.loc[i]) + " |" for i in ct.index] + [""]

    def sv(x):
        return x if isinstance(x, str) else ""
    for c in sorted(d.category.unique()):
        sub = d[(d.category == c) & d.grade_v2.isin(["A_library", "B_multi_tool", "C_one_tool"])].sort_values("height", ascending=False)
        md += [f"## {c} — 이름 근거 A~C ({len(sub)})", "",
               "| feature | m/z | RT | 등급 | 일치 도구 | 출현 | MS-DIAL | MS2Query | SIRIUS | 표시 |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for f, r in sub.head(40).iterrows():
            md.append(f"| {f} | {r.mz:.4f} | {r.rt:.2f} | {r.grade_v2} | {r.tools_agree} | {r.occurrence} | "
                      f"{sv(r.msdial_name)} | {sv(r.q_analog_compound_name) if bool(q_same.get(f, False)) else ''} | "
                      f"{sv(r.s_name)} | {r.name_policy} |")
        md.append("")
    (run / "targets" / f"{sid}_targets_graded_v2.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(info, ensure_ascii=False, indent=1, default=str))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--v2", action="store_true", help="SIRIUS + 출현 근거 판(산출 이름 _v2)")
    ap.add_argument("--sirius-dir", default=None, help="SIRIUS summaries 폴더(기본 runs/<연구>/sirius/summary)")
    args = ap.parse_args()
    cfg = load_study(args.config)
    run, sid = cfg["run_dir"], cfg["study_id"]
    tag = "targets_consensus"

    t = pd.read_csv(run / "targets" / f"{sid}_targets.csv", dtype={"feature_id": str}).set_index("feature_id")
    mg = pd.read_csv(run / "results" / f"{sid}_merged.csv", dtype={"mzmine_id": str}, low_memory=False)
    mg = mg[mg.mzmine_id.notna()].drop_duplicates("mzmine_id").set_index("mzmine_id")[
        ["msdial_name", "msdial_inchikey", "mzmine_name", "msdial_suggestion"]]
    q = pd.read_csv(run / "ms2query" / f"{sid}_ms2query_{tag}_results.csv", dtype={"feature_id": str})
    q = q.groupby("feature_id", sort=False).head(1).set_index("feature_id")[
        ["ms2query_model_prediction", "precursor_mz_difference", "analog_compound_name", "inchikey", "smiles",
         "npc_pathway_results", "npc_superclass_results", "npc_class_results", "cf_class"]].add_prefix("q_")
    m = pd.read_csv(run / "msfinder" / tag / f"{sid}_msfinder_summary.csv", dtype={"feature_id": str}).set_index("feature_id")[
        ["msf_top_name", "msf_top_formula", "msf_top_inchikey", "msf_top_ontology", "msf_top_in_relatives",
         "msf_best_relatives_name", "msf_n_same_formula_top10"]]
    # 2026-11-06: 원결과는 Dropbox 밖 작업 폴더에 생기므로, 끝나면 runs/<연구>/msfinder/<tag>/raw_results/ 로 사본을 둔다.
    # 다른 PC 에서도 등급을 다시 매길 수 있게 두 곳을 모두 보고 가장 최근 것을 쓴다.
    cands = list((MSF_WORK / sid / tag / "out_full").glob("Structure result-*.txt")) +         list((run / "msfinder" / tag / "raw_results").glob("Structure result-*.txt"))
    if not cands:
        raise SystemExit(f"[중단] MS-FINDER 구조 결과 없음: {MSF_WORK / sid / tag / 'out_full'} 또는 {run / 'msfinder' / tag / 'raw_results'}")
    sfile = sorted(cands, key=lambda p: p.stat().st_mtime)[-1]
    sr = pd.read_csv(sfile, sep="\t", dtype=str, keep_default_na=False, quoting=3)
    top10 = {f: {k[:14] for k in g["InChIKey"].head(10) if k} for f, g in sr.groupby("File name", sort=False)}

    d = t.join(mg).join(q).join(m)
    qk = d["q_inchikey"].fillna("").str[:14]
    q_same = (d.q_precursor_mz_difference.abs() <= SAME_MZ_TOL) & (d.q_ms2query_model_prediction >= SCORE_SAME)
    agree = pd.Series([bool(k) and k in top10.get(f, set()) for f, k in zip(d.index, qk)], index=d.index)
    rel = d.msf_top_in_relatives.fillna(False).astype(bool)
    grade = pd.Series("E_unknown", index=d.index)
    grade[d.q_ms2query_model_prediction >= SCORE_CLASS] = "D_class_only"
    grade[q_same | rel] = "C_one_tool"
    grade[q_same & agree] = "B_two_tools"
    grade[d.msdial_name.notna()] = "A_library"
    d["grade"] = grade
    d["ms2query_msfinder_agree"] = agree
    if args.v2:
        return grade_v2(cfg, run, sid, d, q_same, rel, top10, sfile, args.sirius_dir)
    d["best_name"] = d.msdial_name.where(d.msdial_name.notna(),
                     d.q_analog_compound_name.where(q_same, d.msf_best_relatives_name.where(rel)))
    nm = (d.best_name.fillna("") + " " + d.q_analog_compound_name.fillna("")).str.lower()
    d["name_flag"] = nm.apply(lambda s: "library_bias_suspect" if any(w in s for w in SUSPECT_WORDS) else "")
    d.to_csv(run / "targets" / f"{sid}_targets_graded.csv", encoding="utf-8-sig", index_label="feature_id")

    ct = pd.crosstab(d.category, d.grade)
    cls = (d[d.grade.isin(["A_library", "B_two_tools", "C_one_tool", "D_class_only"])]
           .groupby("category").q_npc_pathway_results.agg(lambda s: s.value_counts().head(3).to_dict()))
    info = {"thresholds": {"same_mz_tol_Da": SAME_MZ_TOL, "score_same": SCORE_SAME, "score_class": SCORE_CLASS},
            "grade_by_category": {c: ct.loc[c].to_dict() for c in ct.index},
            "grade_total": d.grade.value_counts().to_dict(),
            "top_npc_pathways_by_category": cls.to_dict(),
            "library_bias_suspect": int((d.name_flag != "").sum()),
            "msfinder_structure_file": str(sfile)}
    (run / "targets" / f"{sid}_targets_graded.json").write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")

    md = [f"# {sid} — 선발 대상 이름 근거 등급", "", "등급은 근거의 겹 수이지 확정이 아니다.", "",
          "## 범주 × 등급", "",
          "| 범주 | " + " | ".join(ct.columns) + " |", "|---" * (len(ct.columns) + 1) + "|"] + \
         [f"| {c} | " + " | ".join(str(v) for v in ct.loc[c]) + " |" for c in ct.index] + [""]
    for c in sorted(d.category.unique()):
        sub = d[(d.category == c) & d.grade.isin(["A_library", "B_two_tools", "C_one_tool"])].sort_values("height", ascending=False)
        md += [f"## {c} — 이름 근거 A~C ({len(sub)})", "",
               "| feature | m/z | RT | 등급 | 이름 | MS2Query 점수 | 계열 | 의심 |", "|---|---|---|---|---|---|---|---|"]
        for f, r in sub.head(40).iterrows():
            md.append(f"| {f} | {r.mz:.4f} | {r.rt:.2f} | {r.grade} | {r.best_name} | {r.q_ms2query_model_prediction} | "
                      f"{r.q_npc_pathway_results} | {r.name_flag} |")
        md.append("")
    (run / "targets" / f"{sid}_targets_graded.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(info, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
