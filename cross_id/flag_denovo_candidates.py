# -*- coding: utf-8 -*-
"""
S6c 신규 구조 후보 칸 — MSNovelist 가 만든 구조를 '시료에서 나올 수 있는가'로 걸러 표시 (2026-11-12, 이전 연구)
=====================================================================================================
무엇을 하나?
  등급 v2 에서 이름을 믿을 수 없는 봉우리(D_class_only · E_unknown = 네 도구 모두 이름 못 붙임)에 한해,
  SIRIUS 6 의 MSNovelist 가 새로 만든 구조 후보를 읽어 네 겹의 걸러내기 표시를 붙인다.
  🔴 등급(도구 일치 겹 수)에는 끼지 않는다. MSNovelist 는 SIRIUS 구조 지문을 입력으로 쓰므로 표를 주면
     SIRIUS 의견을 두 번 세게 된다. 등급 v2 파일은 읽기만 하고 결과는 새 파일에만 쓴다.
  🔴 어느 겹에서 걸려도 행을 지우지 않는다. 판정과 이유를 칸으로 남긴다. '통과'도 표준물질 확인 전까지 가설이다.

네 겹 (판정 규칙은 결과를 보기 전에 고정 — 2026-11-12 승인 설계)
  ① 원소   : 허용 원소(기본 C·H·N·O·S·P) 밖이면 '제외_원소'. 해조류처럼 Br·Cl 이 정상인 연구는
             설정 denovo_guard.allowed_elements 로 넓힌다. 규소(Si)는 실리콘 오염원이라 '오염의심'.
  ② 기록 대조(구조 뼈대 = InChIKey 앞 14자 일치)
             LOTUS 고정판(2026-04-13)에 있으면 '통과_기록'(출현 표지: 종·속·근연·과·다른 생물).
             의약품 기록부(db_frozen/ 고정판: DrugCentral, 선택 ChEMBL 임상)에만 있고 LOTUS·천연물 표지 모두 없으면 '제외_의약품'.
  ③ 가장 가까운 친척(구조 지문 Morgan 반경 2·2048비트, Tanimoto 유사도)
             LOTUS 전체 중 최근접 천연물 유사도 np_nn, 의약품 중 최근접 drug_nn.
             drug_nn ≥ DRUG_MIN 이고 drug_nn ≥ np_nn + DRUG_MARGIN 이면 '검토_합성의심'.
             np_nn ≥ NP_SUPPORT 이면 '통과_골격근거'. 둘 다 아니면 '검토_근거약함'.
             같은 속·근연 속 천연물 중 최근접(rel_nn)과 그 이름도 함께 적는다(판정에는 쓰지 않고 읽는 사람 몫).
             최근접 천연물의 NPClassifier 경로와 SIRIUS 가 그 봉우리에 매긴 경로가 같은지(npc_agree)도 적는다.
             ⚠️ 둘 다 같은 SIRIUS 계산에서 갈라져 나온 정보라 독립 근거로 세지 않는다.
  ④ 실험실 오염 물질: 폴리에틸렌·폴리프로필렌글라이콜 반복, 프탈레이트, 유기인산 트라이에스터,
             장사슬 황산·설폰산·4급 암모늄 계면활성제, 벤젠설폰아마이드 가소제 → '오염의심'
             (사람이 만든 물질이지만 실험 중 실제로 섞여 들어오므로 '시료 유래 불가'와 구분한다.
              범주 출처: Keller BO et al. Anal Chim Acta 2008;627:71-81).

판정 우선순위(위가 이김): 제외_원소 > 제외_의약품 > 오염의심 > 검토_의약품기록 > 통과_기록 > 검토_합성의심 > 통과_골격근거 > 검토_근거약함
  검토_의약품기록 = 의약품 기록이 있고 LOTUS 보고가 먼 분류군(lotus_other)뿐인 구조(LOTUS 의 합성물 검출 보고 대비).
  🔴 LOTUS 에는 합성물의 '생물 검출' 보고도 실려 있다 → 최근접 천연물 기준 목록에서 분류군이 먼 의약품·오염원·허용 밖 원소 구조를 뺀다.

입력
  runs/<sid>/targets/<sid>_targets_graded_v2.csv
  runs/<sid>/sirius/summary/denovo_structure_identifications.tsv   (SIRIUS summaries 가 쓰는 MSNovelist 결과)
  ~/.npdb_cache/lotus_260413_frozen_metadata.csv.gz              (grade_target_annotations.py 와 같은 고정판)
  db_frozen/drugs_*.csv.gz                                       (build_exclusion_db.py, 출처별 가장 최근 판)
산출: runs/<sid>/targets/<sid>_denovo_candidates.{csv,json,md}  (+ 봉우리별 요약 _denovo_by_feature.csv)
쓰는 법:
  py scripts/flag_denovo_candidates.py --config config/<연구>.json                 # 판정
  py scripts/flag_denovo_candidates.py --config config/<연구>.json --write-mgf     # D·E 봉우리만 담은 mgf 쓰기(SIRIUS 입력용)
  --denovo <tsv> 로 다른 위치의 결과를 읽을 수 있다(시험용).
"""
import sys
import re
import json
import argparse
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import utf8_console, load_study

HERE = Path(__file__).resolve().parent
LOTUS_GZ = Path.home() / ".npdb_cache" / "lotus_260413_frozen_metadata.csv.gz"
FROZEN = HERE.parent / "db_frozen"

# 기본 설정값 — 연구 설정의 denovo_guard 로 덮을 수 있다. 결과 json 에 실제 쓴 값을 남긴다.
DEFAULTS = {
    "feature_scope": ["D_class_only", "E_unknown"],   # 네 도구 모두 믿을 이름을 못 붙인 봉우리
    "top_k": 5,                  # 봉우리마다 MSNovelist 상위 몇 개 구조까지 볼지
    "allowed_elements": ["C", "H", "N", "O", "S", "P"],
    "np_support": 0.6,           # 최근접 천연물 유사도가 이 이상이면 '골격 근거 있음'
    "drug_min": 0.5,             # 최근접 약물 유사도가 이 이상이면서
    "drug_margin": 0.1,          #   천연물보다 이만큼 더 가까우면 '합성 의심'
}
CONTAMINANT_SMARTS = {           # 이름 → 부분 구조. Keller 2008 의 오염원 범주를 구조식으로 옮긴 것
    "polyethylene_glycol": "[CH2][CH2]O[CH2][CH2]O[CH2][CH2]O",
    "polypropylene_glycol": "[CH3][CH1]([CH2]O[CH1]([CH3])[CH2]O)O",
    "phthalate_diester": "O=C(O[#6])c1ccccc1C(=O)O[#6]",
    "phthalate_monoester": "O=C([OX2H1,OX1-])c1ccccc1C(=O)O[CH2][CH1,CH2]",   # 가소제 가수분해물(MEHP 등), 2026-11-12 시험에서 추가
    "organophosphate_triester": "P(=O)(O[#6])(O[#6])O[#6]",
    "long_chain_sulfate_sulfonate": "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]~*.S(=O)(=O)[OX2H1,OX1-]",
    "quaternary_ammonium_surfactant": "[NX4+]([CH3])([CH3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]",
    "benzenesulfonamide_plasticizer": "c1ccccc1S(=O)(=O)[NX3][CH2][CH2][CH2][CH3]",
}
PRIORITY = ["제외_원소", "제외_의약품", "오염의심", "검토_의약품기록", "통과_기록", "검토_합성의심", "통과_골격근거", "검토_근거약함"]
CLOSE = ("species", "genus", "relatives", "family")   # 같은 종·속·근연 속·과에서 보고 = 분류군이 가까운 출현 기록


def k14(x):
    return x[:14] if isinstance(x, str) and len(x) >= 14 else ""


def load_drugs():
    """db_frozen/ 의 의약품 기록부를 출처별 가장 최근 판 하나씩 읽어 합친다."""
    files = sorted(FROZEN.glob("drugs_*.csv.gz"))
    if not files:
        raise SystemExit("[중단] 의약품 기록부 고정판 없음 — py scripts/build_exclusion_db.py 를 먼저 돌린다")
    latest = {}
    for f in files:                                   # 이름 = drugs_<출처>_<판>_<날짜> → 출처별 마지막
        latest[f.name.split("_")[1]] = f
    ch = pd.concat([pd.read_csv(f, dtype=str) for f in latest.values()], ignore_index=True)
    return ch, [f.name for f in latest.values()]


def fingerprints(smiles, gen, Chem):
    """SMILES 목록 → (지문 목록, 읽힌 위치). 읽히지 않는 구조는 건너뛴다."""
    fps, keep = [], []
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is not None:
            fps.append(gen.GetFingerprint(m))
            keep.append(i)
    return fps, keep


def nearest(fp, ref_fps, DataStructs):
    """기준 지문 목록 안에서 가장 닮은 것의 (유사도, 위치)."""
    if not ref_fps:
        return 0.0, -1
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    j = max(range(len(sims)), key=sims.__getitem__)
    return sims[j], j


def write_mgf(run, sid, ids):
    """대상 mgf 에서 지정 봉우리만 뽑아 SIRIUS 입력 mgf 를 만든다(FEATURE_ID 줄로 짝을 맞춘다)."""
    src = run / "targets" / f"{sid}_targets.mgf"
    out = run / "targets" / f"{sid}_targets_denovo_scope.mgf"
    keep, block, n = [], [], 0
    for line in src.read_text(encoding="utf-8", errors="replace").splitlines():
        block.append(line)
        if line.strip() == "END IONS":
            fid = next((b.split("=", 1)[1].strip() for b in block if b.upper().startswith("FEATURE_ID=")), None)
            if fid in ids:
                keep.extend(block + [""])
                n += 1
            block = []
    out.write_text("\n".join(keep), encoding="utf-8")
    print(f"저장: {out}  ({n} 봉우리)")


def main():
    utf8_console()
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--denovo", default=None, help="MSNovelist 결과 tsv (기본: sirius/summary/denovo_structure_identifications.tsv)")
    ap.add_argument("--write-mgf", action="store_true")
    args = ap.parse_args()
    cfg = load_study(args.config)
    run, sid = cfg["run_dir"], cfg["study_id"]
    P = {**DEFAULTS, **(cfg.get("denovo_guard") or {})}

    g = pd.read_csv(run / "targets" / f"{sid}_targets_graded_v2.csv", dtype={"feature_id": str}, low_memory=False)
    scope = g[g.grade_v2.isin(P["feature_scope"])].set_index("feature_id")
    if args.write_mgf:
        write_mgf(run, sid, set(scope.index))
        return

    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    dpath = Path(args.denovo) if args.denovo else run / "sirius" / "summary" / "denovo_structure_identifications.tsv"
    dn = pd.read_csv(dpath, sep="\t", dtype=str, keep_default_na=False)
    if dn.empty:
        raise SystemExit(f"[중단] MSNovelist 결과가 비어 있음: {dpath} — SIRIUS 에서 MSNovelist 를 먼저 돌린다")
    dn["rank"] = pd.to_numeric(dn.structurePerIdRank, errors="coerce")
    dn = dn[dn.mappingFeatureId.isin(scope.index) & (dn["rank"] <= P["top_k"])].copy()
    if dn.empty:
        raise SystemExit("[중단] MSNovelist 결과의 mappingFeatureId 가 대상(D·E) 봉우리와 하나도 맞지 않음 — 입력 mgf 확인")

    # ── 기록부 읽기
    lot = pd.read_csv(LOTUS_GZ, dtype=str, usecols=lambda c: c in {
        "structure_inchikey", "structure_smiles_2D", "structure_nameTraditional",
        "structure_taxonomy_npclassifier_01pathway"})
    lot["k14"] = lot.structure_inchikey.str[:14]
    lot = lot.dropna(subset=["k14"]).drop_duplicates("k14").reset_index(drop=True)
    occ = {}                                            # 뼈대 → 가장 가까운 출현 표지(등급 v2 와 같은 순서)
    for tag in ("genus", "relatives", "family"):
        f = run / "db" / f"{sid}_candidates_{tag}.csv"
        if f.exists():
            for k in pd.read_csv(f, dtype=str, usecols=["inchikey"]).inchikey.dropna().str[:14]:
                occ.setdefault(k, tag)
    for k in lot.k14:
        occ.setdefault(k, "lotus_other")
    rel_keys = {k for k, t in occ.items() if t in ("genus", "relatives")}

    ch, drug_files = load_drugs()
    ch["k14"] = ch.inchikey.str[:14]
    ch_by = ch.drop_duplicates("k14").set_index("k14")

    print("지문 계산: LOTUS", len(lot), "· 의약품", len(ch))
    patt = {n: Chem.MolFromSmarts(s) for n, s in CONTAMINANT_SMARTS.items()}
    drug_k14 = set(ch.k14.dropna())
    # 🔴 LOTUS 에는 '어느 생물에서 검출됐다'는 보고로 의약품·가소제 같은 합성물도 실려 있다(2026-11-12 시험:
    #    이부프로펜·시프로플록사신·DEHP·트라이에틸렌글라이콜이 모두 LOTUS 에 있음). 그대로 두면 약물을 닮은 구조가
    #    '천연물을 닮았다'로 통과한다. 그래서 최근접 천연물을 찾는 기준 목록에서는 분류군이 가까운 출현(CLOSE)이 아닌 한
    #    의약품 기록·오염원 구조·허용 밖 원소를 가진 LOTUS 구조를 뺀다.
    lot_fps, lot_keep, n_drop = [], [], {"drug": 0, "contaminant": 0, "element": 0}
    allowed = set(P["allowed_elements"])
    for i, (smi, k) in enumerate(zip(lot.structure_smiles_2D, lot.k14)):
        m0 = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if m0 is None:
            continue
        if occ.get(k) not in CLOSE:
            if k in drug_k14:
                n_drop["drug"] += 1; continue
            if any({a.GetSymbol() for a in m0.GetAtoms()} - allowed):
                n_drop["element"] += 1; continue
            if any(p is not None and m0.HasSubstructMatch(p) for p in patt.values()):
                n_drop["contaminant"] += 1; continue
        lot_fps.append(gen.GetFingerprint(m0))
        lot_keep.append(i)
    lot_rows = lot.iloc[lot_keep].reset_index(drop=True)
    print("최근접 천연물 기준에서 뺀 LOTUS 구조:", n_drop)
    rel_idx = [i for i, k in enumerate(lot_rows.k14) if k in rel_keys]
    rel_fps = [lot_fps[i] for i in rel_idx]
    ch_fps, ch_keep = fingerprints(ch.smiles.tolist(), gen, Chem)
    ch_rows = ch.iloc[ch_keep].reset_index(drop=True)

    out = []
    for d in dn.to_dict("records"):          # 칸 이름에 ':' 가 있어 dict 로 읽는다
        m = Chem.MolFromSmiles(d["smiles"]) if d["smiles"] else None
        row = {"feature_id": d["mappingFeatureId"], "rank": int(d["rank"]), "smiles": d["smiles"],
               "inchikey14": k14(d["InChIkey2D"]), "formula": d["molecularFormula"], "adduct": d["adduct"],
               "model_score": d.get("ModelScore"), "csi_score": d.get("CSI:FingerIDScore")}
        f = scope.loc[d["mappingFeatureId"]]
        row.update({"mz": f.mz, "rt": f.rt, "category": f.category, "grade_v2": f.grade_v2,
                    "s_npc_pathway": f.get("s_npc_pathway", "")})
        reasons = []
        if m is None:
            row.update(status="검토_근거약함", reasons="구조 읽기 실패")
            out.append(row)
            continue
        elems = {a.GetSymbol() for a in m.GetAtoms()} | {"H"}
        bad = sorted(elems - set(P["allowed_elements"]) - {"Si"})
        contam = [n for n, p in patt.items() if p is not None and m.HasSubstructMatch(p)]
        if "Si" in elems:
            contam.append("silicon_siloxane")
        k = row["inchikey14"]
        in_lotus = occ.get(k, "")
        ch_hit = ch_by.loc[k] if k in ch_by.index else None
        fp = gen.GetFingerprint(m)
        np_nn, j = nearest(fp, lot_fps, DataStructs)
        rel_nn, jr = nearest(fp, rel_fps, DataStructs)
        drug_nn, jd = nearest(fp, ch_fps, DataStructs)
        nn = lot_rows.iloc[j] if j >= 0 else None
        nnr = lot_rows.iloc[rel_idx[jr]] if jr >= 0 else None
        nnd = ch_rows.iloc[jd] if jd >= 0 else None
        row.update({
            "elements": "".join(sorted(elems)), "bad_elements": ",".join(bad), "contaminant": ",".join(contam),
            "lotus_occurrence": in_lotus or "none",
            "drug_exact": "" if ch_hit is None else f"{ch_hit.source}:{ch_hit.drug_id}|{ch_hit['name']}",
            "np_nn": np_nn, "np_nn_name": None if nn is None else nn.structure_nameTraditional,
            "np_nn_inchikey14": None if nn is None else nn.k14,
            "np_nn_npc_pathway": None if nn is None else nn.structure_taxonomy_npclassifier_01pathway,
            "rel_nn": rel_nn, "rel_nn_name": None if nnr is None else nnr.structure_nameTraditional,
            "drug_nn": drug_nn, "drug_nn_name": None if nnd is None else f"{nnd['name']} ({nnd.source}:{nnd.drug_id})",
        })
        sp = str(row["s_npc_pathway"] or "")
        npp = str(row["np_nn_npc_pathway"] or "")
        row["npc_agree"] = "" if not sp or not npp or sp == "nan" or npp == "nan" else str(bool(set(sp.split("|")) & set(npp.split("|"))))
        # ── 판정(우선순위 순서로 해당하는 것을 모두 모은 뒤 첫째를 고른다)
        hits = []
        if bad:
            hits.append("제외_원소"); reasons.append(f"허용 밖 원소 {','.join(bad)}")
        if ch_hit is not None and not in_lotus and str(ch_hit.natural_product) != "1":
            hits.append("제외_의약품"); reasons.append(f"의약품 기록만 있음 {ch_hit['name']}")
        if contam:
            hits.append("오염의심"); reasons.append(f"오염원 구조 {','.join(contam)}")
        if in_lotus and in_lotus not in CLOSE and ch_hit is not None:
            hits.append("검토_의약품기록"); reasons.append(f"의약품 기록 + 먼 분류군 LOTUS 보고({in_lotus}) {ch_hit['name']}")
        elif in_lotus:
            hits.append("통과_기록"); reasons.append(f"LOTUS {in_lotus}")
        if drug_nn >= P["drug_min"] and drug_nn >= np_nn + P["drug_margin"]:
            hits.append("검토_합성의심"); reasons.append("천연물보다 약물에 더 닮음")
        if np_nn >= P["np_support"]:
            hits.append("통과_골격근거"); reasons.append("닮은 천연물 있음")
        hits.append("검토_근거약함")
        row["status"] = min(hits, key=PRIORITY.index)
        row["reasons"] = "; ".join(reasons)
        out.append(row)

    res = pd.DataFrame(out).sort_values(["feature_id", "rank"])
    tdir = run / "targets"
    res.to_csv(tdir / f"{sid}_denovo_candidates.csv", index=False, encoding="utf-8-sig")
    # 봉우리별: 통과한 후보 가운데 가장 높은 순위 하나 + 판정 분포
    by = []
    for fid, grp in res.groupby("feature_id"):
        ok = grp[grp.status.str.startswith("통과")]
        best = ok.iloc[0] if len(ok) else None
        by.append({"feature_id": fid, "mz": grp.mz.iloc[0], "rt": grp.rt.iloc[0], "category": grp.category.iloc[0],
                   "grade_v2": grp.grade_v2.iloc[0], "n_candidates": len(grp),
                   "status_counts": json.dumps(grp.status.value_counts().to_dict(), ensure_ascii=False),
                   "best_pass_rank": None if best is None else best["rank"],
                   "best_pass_smiles": None if best is None else best.smiles,
                   "best_pass_status": None if best is None else best.status,
                   "best_pass_np_nn_name": None if best is None else best.np_nn_name})
    byf = pd.DataFrame(by)
    byf.to_csv(tdir / f"{sid}_denovo_by_feature.csv", index=False, encoding="utf-8-sig")
    meta = {"study": sid, "settings": P, "contaminant_smarts": CONTAMINANT_SMARTS, "priority": PRIORITY,
            "inputs": {"graded_v2": str(tdir / f"{sid}_targets_graded_v2.csv"), "denovo": str(dpath),
                       "lotus": LOTUS_GZ.name, "drugs": drug_files},
            "fingerprint": "Morgan radius 2, 2048 bit, Tanimoto (RDKit)",
            "lotus_reference_dropped": n_drop, "lotus_reference_n": len(lot_rows),
            "n_scope_features": int(len(scope)), "n_features_with_denovo": int(res.feature_id.nunique()),
            "n_candidates": int(len(res)), "status_total": res.status.value_counts().to_dict(),
            "features_with_any_pass": int(byf.best_pass_rank.notna().sum())}
    (tdir / f"{sid}_denovo_candidates.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    md = [f"# {sid} — MSNovelist 신규 구조 후보 걸러내기", "",
          "등급 v2 에서 이름을 못 믿는 봉우리(D·E)만 대상. 등급에는 끼지 않는다. '통과'도 표준물질 확인 전까지 가설.", "",
          "| 판정 | 후보 수 |", "|---|---|"] + [f"| {k} | {v} |" for k, v in meta["status_total"].items()] + [
          "", f"통과 후보가 하나라도 있는 봉우리: {meta['features_with_any_pass']} / {meta['n_features_with_denovo']}", ""]
    (tdir / f"{sid}_denovo_candidates.md").write_text("\n".join(md), encoding="utf-8")
    print("저장:", tdir / f"{sid}_denovo_candidates.csv")
    print(json.dumps(meta["status_total"], ensure_ascii=False))


if __name__ == "__main__":
    main()
