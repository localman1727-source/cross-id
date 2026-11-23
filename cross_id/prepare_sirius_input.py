# -*- coding: utf-8 -*-
"""
SIRIUS 입력 고르기 — 공개 사전이 이름을 못 붙인 봉우리 중 믿을 만한 것만 추린다
=================================================================================
'왜' 전부 넣지 않나?
  mzmine 이 내보낸 sirius.mgf 에는 조각무늬가 있는 봉우리가 수만 개다. SIRIUS 는 봉우리마다
  분자식·조각 나무를 풀고 구조 지문을 예측하므로 전부 넣으면 며칠이 걸리고, 잡음·동위원소
  봉우리에 시간을 쓴다. 그래서
    1) 설정의 질량 범위 안
    2) 반복확인용(QC) 측정 세 번 모두에서 잡히고 흔들림이 작다(compare_samples.py 와 같은 기준)
    3) MS-DIAL 이 조각무늬로 이름을 붙이지 못했다(merged 의 msdial_name 이 비어 있다)
    4) 위 조건을 만족하는 것 중 신호가 큰 순서로 --top 개
  를 고른다. 4) 는 계산 시간 때문에 두는 상한이지 생물학적 기준이 아니다 — 고른 기준은
  selection.json 에 남긴다.

입력  C:/mzmine_out/<연구>/<연구>_<극성>_sirius.mgf   (make_mzmine_batches.py 의 OUT_BIG)
      runs/<연구>/mzmine_out/*_full_feature_table.csv
      runs/<연구>/results/<연구>_merged.csv            (merge_annotations.py, MS-DIAL 결과 포함판)
산출  runs/<연구>/sirius/<연구>_sirius_input.mgf · _selection.csv · _selection.json
실행  py scripts/prepare_sirius_input.py --config config/<연구>.json [--top 2000]
"""
import sys
import csv
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console
from compare_samples import QC_RSD_MAX

utf8_console()
OUT_BIG_ROOT = Path(r"C:\mzmine_out")        # make_mzmine_batches.py 와 같은 자리


def read_mgf_blocks(path):
    """FEATURE_ID → [블록 문자열…], 블록별 MSLEVEL."""
    blocks, cur = {}, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            cur.append(line)
            if line.startswith("END IONS"):
                txt = "".join(cur)
                cur = []
                fid = next((l.split("=", 1)[1].strip() for l in txt.splitlines()
                            if l.startswith("FEATURE_ID=")), None)
                lvl = next((l.split("=", 1)[1].strip() for l in txt.splitlines()
                            if l.startswith("MSLEVEL=")), None)
                if fid is not None:
                    blocks.setdefault(fid, []).append((lvl, txt))
    return blocks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--top", type=int, default=2000, help="신호가 큰 순서로 이만큼만 (계산 시간 상한)")
    args = ap.parse_args()

    cfg = load_study(args.config)
    run = cfg["run_dir"]
    sid = cfg["study_id"]
    pol = cfg["polarity"].upper()
    pr = cfg["processing"]

    mgf = OUT_BIG_ROOT / sid / f"{sid}_{pol}_sirius.mgf"
    merged_csv = run / "results" / f"{sid}_merged.csv"
    feat_csv = next((run / "mzmine_out").glob("*_full_feature_table.csv"))
    for p in (mgf, merged_csv):
        if not p.exists():
            raise SystemExit(f"[중단] 없음: {p}")

    # 2) QC 세 번 모두 검출 + 흔들림 상한
    ft = pd.read_csv(feat_csv, low_memory=False)
    qc_cols = [c for c in ft.columns if c.startswith("datafile:") and c.endswith(":height") and "_QC_" in c]
    qc = ft[qc_cols].apply(pd.to_numeric, errors="coerce")
    qc_all = qc.notna().all(axis=1) & (qc > 0).all(axis=1)
    qc_rsd = qc.std(axis=1, ddof=1) / qc.mean(axis=1)
    ft["qc_ok"] = qc_all & (qc_rsd <= QC_RSD_MAX)
    ft["qc_rsd"] = qc_rsd
    ft["id"] = ft["id"].astype(str)

    # 3) MS-DIAL 이 조각무늬로 이름 붙인 봉우리 제외
    mg = pd.read_csv(merged_csv, low_memory=False, dtype={"mzmine_id": str})
    if "mzmine_id" not in mg.columns:
        raise SystemExit("[중단] merged 에 mzmine_id 가 없습니다 — merge_annotations.py 를 다시 돌리세요")
    if not (mg["found_by"] != "mzmine 만").any():
        raise SystemExit("[중단] merged 에 MS-DIAL 결과가 없습니다 — MS-DIAL 실행 뒤 merge 를 다시 돌리세요")
    named_by_msdial = set(mg.loc[mg["msdial_name"].notna() & mg["mzmine_id"].notna(), "mzmine_id"])

    blocks = read_mgf_blocks(mgf)
    rows = []
    for fid, bl in blocks.items():
        if not any(l == "2" for l, _ in bl):
            continue
        rows.append({"id": fid})
    cand = pd.DataFrame(rows).merge(ft[["id", "mz", "rt", "height", "qc_ok", "qc_rsd"]], on="id", how="left")
    n_ms2 = len(cand)
    cand = cand[(cand["mz"] >= pr["mz_min"]) & (cand["mz"] <= pr["mz_max"])]
    n_range = len(cand)
    cand = cand[cand["qc_ok"] == True]
    n_qc = len(cand)
    cand = cand[~cand["id"].isin(named_by_msdial)]
    n_unnamed = len(cand)
    sel = cand.sort_values("height", ascending=False).head(args.top)

    out = run / "sirius"
    out.mkdir(parents=True, exist_ok=True)
    with (out / f"{sid}_sirius_input.mgf").open("w", encoding="utf-8", newline="\n") as fh:
        for fid in sel["id"]:
            for _, txt in blocks[fid]:
                fh.write(txt)
    sel.to_csv(out / f"{sid}_sirius_selection.csv", index=False, encoding="utf-8-sig")
    info = {"mgf_source": str(mgf), "criteria": {"mz_min": pr["mz_min"], "mz_max": pr["mz_max"],
            "qc_detected_all": True, "qc_rsd_max": QC_RSD_MAX,
            "exclude_msdial_spectrum_matched": True, "rank": "max MS1 height (mzmine 'height')",
            "top": args.top},
            "counts": {"features_with_ms2": n_ms2, "in_mass_range": n_range, "qc_stable": n_qc,
                       "not_named_by_msdial": n_unnamed, "selected": len(sel)}}
    (out / f"{sid}_sirius_selection.json").write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(info["counts"], ensure_ascii=False))
    print(f"[생성] {out}")


if __name__ == "__main__":
    main()
