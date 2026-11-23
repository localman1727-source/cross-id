# -*- coding: utf-8 -*-
"""
S6b 이름의 천연물 출처 표시 — 등급 v2 에 붙은 이름이 '천연물 목록에 실린 구조'인가?
========================================================================
'왜' 필요한가? (2026-11-09 이전 연구 결과)
  SIRIUS 를 BIO 묶음만 찾도록 했는데도 BIO 묶음에 의약품·독성 목록(DSSTox·PubChem 의약품 표지)이 섞여 있어,
  두 도구가 일치한 B 등급 이름 가운데 천연물로 보기 어려운 구조(프로테아제 저해제·항생제류)가 상당수 나왔다.
  등급(도구 일치 겹 수)은 그대로 두고, 그 이름의 구조가 천연물 목록에 실렸는지를 '나란히' 표시한다.

천연물 근거로 세는 것 (하나라도 있으면 np_support=True)
  ① LOTUS 출현 기록 — 등급 v2 의 occurrence 가 none 이 아님 (어느 생물에서든 보고됨)
  ② SIRIUS 1순위 구조가 일치 구조일 때, 그 links 칸의 천연물 목록 NP_DBS_SIRIUS
  ③ MS-FINDER 상위 10 안 일치 구조의 Databases 칸의 천연물 목록 NP_DBS_MSF
  np_drug_conflict = 천연물 목록에도 있으나 의약품·독성 표지가 붙고 LOTUS 기록이 없는 구조(반합성 항생제 등) → 사람 검토 우선.
  SuperNatural 은 반합성·유사체를 포함하므로 혼자서는 세지 않고 칸(supernatural_only)으로만 남긴다.
  MS2Query 이름은 사전 출처 칸이 없어 ①로만 판단된다.

산출: targets/<sid>_targets_graded_v2_npsource.{csv,json,md}  — 등급 v2 파일은 건드리지 않는다.
쓰는 법: py scripts/flag_np_source.py --config config/<연구>.json
"""
import sys
import re
import json
import argparse
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import utf8_console, load_study

MSF_WORK = Path(r"C:\MSFINDER\work")       # 윈도우 전용 도구 작업 폴더(절대경로 예외, grade_target_annotations.py 와 같음)
NP_DBS_SIRIUS = {"COCONUT", "LOTUS", "KNAPSACK", "PLANTCYC", "TeroMol", "MiMeDB"}
NP_DBS_MSF = {"COCONUT", "UNPD", "KNApSAcK", "KNAPSACK", "NANPDB", "PlantCyc", "PLANTCYC", "LOTUS", "STOFF"}
DRUG_TAGS = {"PUBCHEMANNOTATIONDRUG", "DSSTox", "PUBCHEMANNOTATIONSAFETYANDTOXIC"}
NAMED = ["A_library", "B_multi_tool", "C_one_tool"]


def sirius_dbs(links):
    """links 칸 'DB:(id);DB2:(id)' 와 'PUBCHEMANNOTATIONx' 표지에서 목록 이름 집합."""
    if not isinstance(links, str):
        return set()
    return {p.split(":")[0].strip() for p in links.split(";") if p.strip()}


def main():
    utf8_console()
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_study(args.config)
    run, sid = cfg["run_dir"], cfg["study_id"]
    tag = "targets_consensus"

    g = pd.read_csv(run / "targets" / f"{sid}_targets_graded_v2.csv", dtype={"feature_id": str}, low_memory=False)
    s = pd.read_csv(run / "sirius" / "summary" / "structure_identifications.tsv", sep="\t", dtype={"mappingFeatureId": str})
    s = s.sort_values("structurePerIdRank").drop_duplicates("mappingFeatureId").set_index("mappingFeatureId")
    s_db = s["links"].map(sirius_dbs)

    cands = list((MSF_WORK / sid / tag / "out_full").glob("Structure result-*.txt")) + \
        list((run / "msfinder" / tag / "raw_results").glob("Structure result-*.txt"))
    if not cands:
        raise SystemExit("[중단] MS-FINDER 구조 결과 없음")
    sfile = sorted(cands, key=lambda p: p.stat().st_mtime)[-1]
    sr = pd.read_csv(sfile, sep="\t", dtype=str, keep_default_na=False, quoting=3)
    msf = {}   # (feature, inchikey14) -> 목록 이름 집합 (상위 10 안)
    for f, grp in sr.groupby("File name", sort=False):
        for ik, dbs in zip(grp["InChIKey"].head(10), grp["Databases"].head(10)):
            if ik:
                msf.setdefault((f, ik[:14]), set()).update(re.findall(r"([A-Za-z]+)=", dbs))

    rows = []
    for r in g.itertuples(index=False):
        key = r.agree_inchikey14 if isinstance(r.agree_inchikey14, str) else ""
        tools = r.tools_agree if isinstance(r.tools_agree, str) else ""
        sd = s_db.get(r.feature_id, set()) if (key and r.s_inchikey14 == key) else set()
        md = msf.get((r.feature_id, key), set()) if key else set()
        occ = r.occurrence if isinstance(r.occurrence, str) else "none"
        s_np = sorted(sd & NP_DBS_SIRIUS)
        m_np = sorted(md & NP_DBS_MSF)
        np_support = (occ != "none") or bool(s_np) or bool(m_np)
        rows.append({"feature_id": r.feature_id, "np_key_inchikey14": key,
                     "np_lotus_occurrence": occ,
                     "np_sirius_dbs": ";".join(s_np), "np_msfinder_dbs": ";".join(m_np),
                     "sirius_drug_tox_tags": ";".join(sorted(sd & DRUG_TAGS)),
                     "supernatural_only": (not np_support) and ("SUPERNATURAL" in sd),
                     # 천연물 목록에도 있지만 의약품·독성 표지도 붙고 LOTUS 출현 기록은 없는 구조(반합성 항생제 등) — 사람 검토 우선
                     "np_drug_conflict": np_support and occ == "none" and bool(sd & DRUG_TAGS),
                     "lipidmaps": "LIPIDMAPS" in sd,
                     "np_support": np_support})
    f = pd.DataFrame(rows)
    out = g.merge(f, on="feature_id", how="left")
    named = out.grade_v2.isin(NAMED)
    # 이름 표시 정책: 이름 근거 등급(A~C)이라도 천연물 근거가 없으면 계열만 내건다(등급 자체는 그대로).
    out["name_policy_np"] = out["name_policy"].where(~named | out.np_support, "class_only_non_np")
    out.to_csv(run / "targets" / f"{sid}_targets_graded_v2_npsource.csv", encoding="utf-8-sig", index=False)

    ct = pd.crosstab(out.grade_v2, out.np_support)
    ct_tools = pd.crosstab(out[named].tools_agree, out[named].np_support)
    info = {"study": sid, "source_grade_file": f"{sid}_targets_graded_v2.csv", "msfinder_structure_file": str(sfile),
            "np_dbs_sirius": sorted(NP_DBS_SIRIUS), "np_dbs_msfinder": sorted(NP_DBS_MSF),
            "drug_tox_tags": sorted(DRUG_TAGS), "supernatural_counted": False,
            "grade_x_np_support": {k: {str(c): int(v) for c, v in ct.loc[k].items()} for k in ct.index},
            "named_tools_x_np_support": {k: {str(c): int(v) for c, v in ct_tools.loc[k].items()} for k in ct_tools.index},
            "named_np_support_by_source": {
                "lotus": int((named & (out.np_lotus_occurrence != "none")).sum()),
                "sirius_np_db": int((named & (out.np_sirius_dbs.fillna("") != "")).sum()),
                "msfinder_np_db": int((named & (out.np_msfinder_dbs.fillna("") != "")).sum())},
            "named_np_drug_conflict": int((named & out.np_drug_conflict.fillna(False)).sum()),
            "named_non_np_lipidmaps": int((named & ~out.np_support & out.lipidmaps.fillna(False)).sum()),
            "named_supernatural_only": int((named & out.supernatural_only.fillna(False)).sum()),
            "named_with_drug_tox_tag_and_no_np": int((named & ~out.np_support & (out.sirius_drug_tox_tags.fillna("") != "")).sum()),
            "name_policy_np_counts": out.name_policy_np.value_counts().to_dict()}
    (run / "targets" / f"{sid}_targets_graded_v2_npsource.json").write_text(
        json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")

    md = [f"# {sid} — 이름의 천연물 출처 표시 (등급 v2 위에 나란히)", "",
          "천연물 근거 = LOTUS 출현 기록 · SIRIUS 일치 구조의 천연물 목록 · MS-FINDER 일치 구조의 천연물 목록 중 하나. SuperNatural 단독은 세지 않음.", "",
          "| 등급 | 천연물 근거 있음 | 없음 |", "|---|---|---|"]
    for k in ct.index:
        md.append(f"| {k} | {int(ct.loc[k].get(True, 0))} | {int(ct.loc[k].get(False, 0))} |")
    md += ["", "## 이름 근거 등급(A~C) 중 천연물 근거 있는 것", "",
           "| feature | m/z | RT | 범주 | 등급 | 일치 도구 | 출현 | 의약품 표지 겹침 | SIRIUS 목록 | MS-FINDER 목록 | SIRIUS 이름 | MS-FINDER 이름 | MS2Query 이름 |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in out[named & out.np_support].sort_values(["grade_v2", "category"]).itertuples():
        md.append(f"| {r.feature_id} | {r.mz:.4f} | {r.rt:.2f} | {r.category} | {r.grade_v2} | {r.tools_agree} | {r.np_lotus_occurrence} | {"검토" if r.np_drug_conflict else ""} | "
                  f"{r.np_sirius_dbs if isinstance(r.np_sirius_dbs, str) else ''} | {r.np_msfinder_dbs if isinstance(r.np_msfinder_dbs, str) else ''} | "
                  f"{r.s_name if isinstance(r.s_name, str) else ''} | {r.msf_top_name if isinstance(r.msf_top_name, str) else ''} | "
                  f"{r.q_analog_compound_name if isinstance(r.q_analog_compound_name, str) else ''} |")
    (run / "targets" / f"{sid}_targets_graded_v2_npsource.md").write_text("\n".join(md), encoding="utf-8")
    print(f"[생성] {run / 'targets' / (sid + '_targets_graded_v2_npsource.csv')} (값은 json 에)")


if __name__ == "__main__":
    main()
