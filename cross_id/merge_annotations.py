# -*- coding: utf-8 -*-
"""
⑥ 두 도구 결과 맞대보기
=========================
'왜' 이름이 아니라 질량과 머무는 시간으로 맞추나?
  같은 성분이라도 두 도구가 붙인 이름이 다를 수 있고, 이름이 같아도 다른 것을 가리킬 수 있다.
  기계가 실제로 잰 값은 질량과 머무는 시간이므로 그것으로 짝을 짓고, 이름은 그다음에 본다.

'왜' 등급을 나누나?
  둘 다 같은 이름을 지목한 성분과 한쪽만 지목한 성분을 한 표에 섞으면, 읽는 사람은 모두
  같은 무게로 받아들인다. 등급을 나눠 두면 어디까지 믿고 쓸지 스스로 고를 수 있다.

등급 (2026-10-08 개정 — 이름 글자가 아니라 구조 지문 InChIKey 로 판정, resolve_structures.py)
  A  두 도구가 같은 봉우리에 같은 구조(또는 같은 뼈대)를 지목   → 두 도구가 일치
  B  두 도구가 같은 봉우리에 다른 구조를 지목                   → 확인 대상
       name_relation 칸: same_formula_different_structure(이성질체 후보끼리)
                         different_formula(이온 형태를 다르게 읽었을 가능성이 크다)
  C  한 도구만 이름을 붙였다                               → 보류
  D  둘 다 찾았지만 이름이 없다                            → 미동정 (많은 것이 정상이다)

⚠️ 이 표의 이름은 확정이 아니다. 같은 분자식의 다른 구조는 원리적으로 구별되지 않는다.

'왜' MS-DIAL 이름은 조각무늬가 맞은 것만 이름으로 치나? (2026-10-08)
  MS-DIAL 은 질량만 맞고 조각무늬 점수가 낮은 후보에도 이름을 적는다. 5.5.260820 은 그 앞에
  "low score:" 를 붙였지만 5.5.261008 부터는 떼고 근거 종류(opt_global_evidence_source)를 따로
  적는다. 이름 칸만 보면 낮은 후보가 이름으로 둔갑해 A·B 등급이 부풀므로,
  opt_global_spectrum_matched=True 인 것만 이름(msdial_name)으로 쓰고 나머지는
  msdial_suggestion 에 남긴다.
  mzmine 쪽 이름은 후보 목록과 질량만 맞춘 것(simple_compound_db_annotation)이다. 한 봉우리에
  후보가 여럿이면 첫 줄을 이름으로 쓰고 후보 수를 mzmine_n_candidates 에 남긴다.

입력
  runs/<study>/mzmine_out/*_annotations.csv          mzmine 동정 결과
  runs/<study>/mzmine_out/*_full_feature_table.csv   mzmine 성분표 (질량·시간)
  runs/<study>/msdial_out/*.mzTab                    MS-DIAL 정렬 결과
산출물
  runs/<study>/results/<study>_merged.csv  ·  <study>_merged.md
실행: py scripts/merge_annotations.py --config config/<연구>.json [--mz-ppm 10 --rt-min 0.2]
"""
import sys
import csv
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console
import resolve_structures as RS

utf8_console()


# ── MS-DIAL mzTab 읽기 ─────────────────────────────────────────────────────
def read_mztab(path: Path) -> list:
    """SMF(잰 값: 질량·시간)와 SML(붙인 이름)을 이어 붙여 한 줄씩 돌려준다.

    '왜' 두 구역을 이어야 하나? mzTab 은 잰 값과 붙인 이름을 다른 구역에 적는다.
    한쪽만 읽으면 이름 없는 목록이거나 위치 없는 이름 목록이 된다.
    """
    smf_head, sml_head = [], []
    smf, sml = {}, []
    for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
        cell = line.split("\t")
        tag = cell[0]
        if tag == "SFH":
            smf_head = cell
        elif tag == "SMF" and smf_head:
            row = dict(zip(smf_head, cell))
            smf[row.get("SMF_ID")] = row
        elif tag == "SMH":
            sml_head = cell
        elif tag == "SML" and sml_head:
            sml.append(dict(zip(sml_head, cell)))

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    out = []
    for r in sml:
        ref = (r.get("SMF_ID_REFS") or "").split("|")[0]
        f = smf.get(ref, {})
        rt_s = num(f.get("retention_time_in_seconds"))
        name = r.get("chemical_name")
        name = None if name in (None, "", "null") else name
        matched = str(r.get("opt_global_spectrum_matched")).strip().lower() == "true"
        out.append({
            "mz": num(f.get("exp_mass_to_charge")),
            "rt": rt_s / 60.0 if rt_s is not None else None,
            "name": name if matched else None,
            "suggestion": None if matched else name,
            "evidence": r.get("opt_global_evidence_source"),
            "formula": (None if r.get("chemical_formula") in (None, "", "null") else r.get("chemical_formula")) if matched else None,
            "adduct": r.get("adduct_ions"),
            "spectrum_matched": r.get("opt_global_spectrum_matched"),
        })
    return out


# ── mzmine 결과 읽기 ───────────────────────────────────────────────────────
def read_mzmine(run_dir: Path) -> list:
    """성분표와 동정표를 합쳐 한 줄씩 돌려준다. 칸 이름이 판마다 조금씩 달라 후보를 넓게 본다."""
    out_dir = run_dir / "mzmine_out"
    feats = sorted(out_dir.glob("*_full_feature_table.csv"))
    annos = sorted(out_dir.glob("*_annotations.csv"))
    if not feats:
        return []

    def pick(row, *names):
        for n in names:
            for k in row:
                if k.lower() == n.lower() or k.lower().endswith(":" + n.lower()):
                    v = row[k]
                    if v not in (None, ""):
                        return v
        return None

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    rows = []
    with feats[0].open(encoding='utf-8-sig') as fh:
        for r in csv.DictReader(fh):
            rows.append({
                "id": pick(r, "id", "row ID"),
                "mz": num(pick(r, "mz", "row m/z")),
                "rt": num(pick(r, "rt", "row retention time")),
                "name": None,
                "formula": None,
                "n_candidates": 0,
            })

    if annos:
        byid, ncand = {}, {}
        with annos[0].open(encoding='utf-8-sig') as fh:
            for a in csv.DictReader(fh):
                key = pick(a, "id", "row ID")
                if key is None:
                    continue
                ncand[key] = ncand.get(key, 0) + 1
                if key not in byid:
                    byid[key] = a
        for r in rows:
            a = byid.get(r["id"])
            if a:
                r["name"] = pick(a, "compound_name", "compound name", "name")
                r["formula"] = pick(a, "formula", "molecular_formula", "mol_formula")
                r["n_candidates"] = ncand.get(r["id"], 0)
    return rows


def merge(mz_rows, md_rows, mz_ppm, rt_tol):
    """질량과 머무는 시간이 모두 창 안에 들면 같은 성분으로 본다."""
    used = set()
    merged = []

    def same_name(a, b):
        if not a or not b:
            return False
        norm = lambda s: "".join(ch for ch in s.lower() if ch.isalnum())
        return norm(a) == norm(b)

    for m in mz_rows:
        best, best_d = None, None
        if m["mz"] is not None:
            tol = m["mz"] * mz_ppm / 1e6
            for i, d in enumerate(md_rows):
                if i in used or d["mz"] is None:
                    continue
                if abs(d["mz"] - m["mz"]) > tol:
                    continue
                if m["rt"] is not None and d["rt"] is not None and abs(d["rt"] - m["rt"]) > rt_tol:
                    continue
                dist = abs(d["mz"] - m["mz"])
                if best_d is None or dist < best_d:
                    best, best_d, best_i = d, dist, i
        if best is not None:
            used.add(best_i)
            if m["name"] and best["name"]:
                grade = "A" if same_name(m["name"], best["name"]) else "B"
            elif m["name"] or best["name"]:
                grade = "C"
            else:
                grade = "D"
            merged.append({"grade": grade, "mz": m["mz"], "rt": m["rt"],
                           "mzmine_name": m["name"], "msdial_name": best["name"],
                           "mzmine_formula": m["formula"], "msdial_formula": best["formula"],
                           "msdial_adduct": best["adduct"], "found_by": "둘 다", "mzmine_id": m["id"],
                           "mzmine_n_candidates": m["n_candidates"],
                           "msdial_suggestion": best["suggestion"], "msdial_evidence": best["evidence"]})
        else:
            merged.append({"grade": "C" if m["name"] else "D", "mz": m["mz"], "rt": m["rt"],
                           "mzmine_name": m["name"], "msdial_name": None,
                           "mzmine_formula": m["formula"], "msdial_formula": None,
                           "msdial_adduct": None, "found_by": "mzmine 만", "mzmine_id": m["id"],
                           "mzmine_n_candidates": m["n_candidates"],
                           "msdial_suggestion": None, "msdial_evidence": None})

    for i, d in enumerate(md_rows):
        if i in used:
            continue
        merged.append({"grade": "C" if d["name"] else "D", "mz": d["mz"], "rt": d["rt"],
                       "mzmine_name": None, "msdial_name": d["name"],
                       "mzmine_formula": None, "msdial_formula": d["formula"],
                       "msdial_adduct": d["adduct"], "found_by": "MS-DIAL 만", "mzmine_id": None,
                       "mzmine_n_candidates": 0,
                       "msdial_suggestion": d["suggestion"], "msdial_evidence": d["evidence"]})
    return merged


def annotate_structures(merged, md_idx, mz_idx, pubchem, rt_tol):
    """줄마다 두 도구 이름의 지문을 붙이고, 지문 관계로 A/B 를 다시 정한다.
    같은 지문이 여러 시간대에 나오면 표시한다 — 한 물질은 한 시간대에만 나오므로,
    여러 시간대면 배당체 등이 기기 안에서 떨어져 같은 조각을 낸 것일 수 있다."""
    for r in merged:
        r["mzmine_inchikey"], r["mzmine_key_source"] = RS.ikey_for(r["mzmine_name"], mz_idx, pubchem)
        r["msdial_inchikey"], r["msdial_key_source"] = RS.ikey_for(r["msdial_name"], md_idx, pubchem)
        if r["mzmine_name"] and r["msdial_name"]:
            rel = RS.relation(r["mzmine_inchikey"], r["msdial_inchikey"],
                              r["mzmine_formula"], r["msdial_formula"],
                              r["mzmine_name"], r["msdial_name"])
            r["name_relation"] = rel
            r["grade"] = "A" if rel in ("same_structure", "same_skeleton_stereo_differs", "same_name_no_key") else "B"
        else:
            r["name_relation"] = None
    # 같은 뼈대가 몇 개의 시간대에 나오나 (rt_tol 보다 떨어진 무리 수)
    rts = {}
    for r in merged:
        for k in (r["msdial_inchikey"], r["mzmine_inchikey"] if r["grade"] == "A" else None):
            if k and r["rt"] is not None:
                rts.setdefault(k[:14], set()).add(r["rt"])
    def n_groups(v):
        v = sorted(v); g = 1
        for a, b in zip(v, v[1:]):
            if b - a > rt_tol:
                g += 1
        return g
    ngrp = {k: n_groups(v) for k, v in rts.items()}
    for r in merged:
        k = r["msdial_inchikey"] or (r["mzmine_inchikey"] if r["grade"] == "A" else None)
        r["same_id_rt_groups"] = ngrp.get(k[:14]) if k else None
        key = r["msdial_inchikey"] or r["mzmine_inchikey"]
        r["pubchem_cid"] = r["pubchem_title"] = None
        if r["grade"] in "AB" and key and pubchem is not None:
            t = pubchem.title_from_ikey(key)
            if t:
                r["pubchem_cid"], r["pubchem_title"] = t.get("cid"), t.get("title")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--mz-ppm", type=float, default=10.0, help="질량이 이만큼 안에서 맞으면 같은 것으로 본다")
    ap.add_argument("--rt-min", type=float, default=0.2, help="머무는 시간이 이만큼(분) 안에서 맞으면 같은 것으로 본다")
    ap.add_argument("--offline", action="store_true", help="PubChem 을 조회하지 않는다(캐시만 쓴다)")
    ap.add_argument("--mztab", default=None, help="다른 mzTab 을 쓸 때(시험 실행 등)")
    ap.add_argument("--out-dir", default=None, help="산출물 폴더를 바꿀 때(기본 runs/<연구>/results)")
    args = ap.parse_args()

    cfg = load_study(args.config)
    run = cfg["run_dir"]

    mz_rows = read_mzmine(run)
    md_files = ([Path(args.mztab)] if args.mztab else
                sorted((run / "msdial_out").glob("*.mzTab")) if (run / "msdial_out").exists() else [])
    md_rows = read_mztab(md_files[0]) if md_files else []

    if not mz_rows and not md_rows:
        raise SystemExit("[중단] 두 도구의 결과가 모두 없습니다 — 먼저 run_stage1.py 로 분석을 돌리세요")
    if not mz_rows:
        print("  ⚠ mzmine 결과가 없습니다 — MS-DIAL 쪽만 정리합니다")
    if not md_rows:
        print("  ⚠ MS-DIAL 결과가 없습니다 — mzmine 쪽만 정리합니다")

    merged = merge(mz_rows, md_rows, args.mz_ppm, args.rt_min)
    msp = cfg["processing"].get("msdial_msp")
    md_idx = RS.msp_index(msp) if msp and Path(msp).exists() else {}
    mz_idx = RS.candidate_index(cfg["db_dir"])
    pubchem = RS.PubChem(offline=args.offline)
    annotate_structures(merged, md_idx, mz_idx, pubchem, args.rt_min)
    pubchem.save()
    print(f"  PubChem 조회 {pubchem.n_calls} 회 (나머지는 캐시)")
    merged.sort(key=lambda r: ("ABCD".index(r["grade"]), -(r["mz"] or 0)))

    res = Path(args.out_dir) if args.out_dir else run / "results"
    res.mkdir(parents=True, exist_ok=True)
    out_csv = res / f"{cfg['study_id']}_merged.csv"
    with out_csv.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(merged[0].keys()))
        w.writeheader(); w.writerows(merged)

    counts = {g: sum(1 for r in merged if r["grade"] == g) for g in "ABCD"}
    lines = [
        f"# {cfg['study_id']} — 두 도구 결과 맞대보기",
        "",
        f"맞춘 기준: 질량 {args.mz_ppm} ppm · 머무는 시간 {args.rt_min} 분 이내.",
        f"mzmine {len(mz_rows)} 줄 · MS-DIAL {len(md_rows)} 줄 → 합친 표 {len(merged)} 줄.",
        "",
        "| 등급 | 뜻 | 개수 |",
        "|---|---|---|",
        f"| A | 같은 봉우리에 두 도구가 같은 구조(지문)를 지목 | {counts['A']} |",
        f"| B | 같은 봉우리에 다른 구조를 지목 (확인 대상) | {counts['B']} |",
        f"| C | 한 도구만 이름을 붙였다 (보류) | {counts['C']} |",
        f"| D | 둘 다 찾았지만 이름이 없다 (미동정) | {counts['D']} |",
        "",
        f"MS-DIAL 이름은 조각무늬가 사전과 맞은 것({sum(1 for d in md_rows if d['name'])} 줄)만 셌다. 낮은 점수 후보 "
        f"{sum(1 for d in md_rows if d['suggestion'])} 줄은 msdial_suggestion 칸에만 있다.",
        "mzmine 이름은 후보 목록과 질량만 맞춘 것이다(조각무늬 대조 없음).",
        "A/B 는 이름 글자가 아니라 구조 지문(InChIKey)으로 판정했다. 지문은 두 도구의 원본 사전에서 붙였고,",
        "없을 때만 PubChem 에서 찾았다(mzmine_key_source / msdial_key_source 칸).",
        f"같은 지문이 두 시간대 이상에 나오는 줄(same_id_rt_groups≥2): "
        f"{sum(1 for r in merged if (r.get('same_id_rt_groups') or 0) >= 2)} — 배당체 등이 기기 안에서 떨어진 조각일 수 있다.",
        "",
        "이 표의 이름은 확정이 아니다. 같은 분자식의 다른 구조는 원리적으로 구별되지 않으므로,",
        "확정하려면 표준품을 같은 조건에서 돌려 머무는 시간까지 맞춰야 한다.",
    ]
    if counts["B"]:
        lines += ["", "## 다른 구조를 지목한 성분", "", "| 질량 | 시간(분) | mzmine | MS-DIAL | 관계 |", "|---|---|---|---|---|"]
        for r in [x for x in merged if x["grade"] == "B"][:30]:
            lines.append(f"| {r['mz']} | {r['rt']} | {r['mzmine_name']} | {r['msdial_name']} | {r['name_relation']} |")

    out_md = res / f"{cfg['study_id']}_merged.md"
    out_md.write_text("\n".join(lines) + "\n", encoding='utf-8')

    print(f"[생성] {out_csv}")
    print(f"[생성] {out_md}")
    print(f"  A {counts['A']} · B {counts['B']} · C {counts['C']} · D {counts['D']}")


if __name__ == "__main__":
    main()
