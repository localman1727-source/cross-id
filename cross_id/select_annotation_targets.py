# -*- coding: utf-8 -*-
"""
이름 붙이기 대상 고르기 — 선발 범주(compare_samples.py 5단계) 기준
=================================================================
'왜' 바꾸나? (2026-10-08 결정)
  처음 대상(prepare_sirius_input.py)은 'MS-DIAL 이 이름 못 붙임 + QC 안정 + 신호 큰 순 2000'이었다.
  이 연구가 찾는 것은 산지 사이에서 다르거나(산지 특이·우세) 세 산지에 공통인 물질이므로, 대상도
  그 범주로 고른다. 신호 크기 순은 계산 시간 상한일 뿐이었다.
    1) <산지>_only · <산지>_dominant 는 조각무늬(MS2)가 있으면 전부
    2) common_similar · common_present 중 QC 안정(qc_rsd ≤ QC_RSD_MAX)인 것을 신호 큰 순 --common-top 개
  MS-DIAL 이 이미 이름 붙인 봉우리도 빼지 않는다 — 두 도구 결과가 같은 쪽을 가리키는지 보는 근거가 된다.
  질량 범위(설정 processing.mz_min~mz_max)는 prepare_sirius_input.py 와 같게 둔다.

입력  runs/<연구>/compare/<연구>_<극성>_selection_categories.csv
      C:/mzmine_out/<연구>/<연구>_<극성>_sirius.mgf   (런별 MS1·MS2 전부)
산출  runs/<연구>/targets/<연구>_targets.csv · _targets.mgf · _targets.json
실행  py scripts/select_annotation_targets.py --config config/<연구>.json [--common-top 1000]
"""
import sys
import json
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console
from compare_samples import QC_RSD_MAX

utf8_console()
OUT_BIG_ROOT = Path(r"C:\mzmine_out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--common-top", type=int, default=1000, help="공통 범주에서 신호 큰 순 이만큼 (계산 시간 상한)")
    args = ap.parse_args()
    cfg = load_study(args.config)
    run, sid, pol, pr = cfg["run_dir"], cfg["study_id"], cfg["polarity"].upper(), cfg["processing"]

    cats = pd.read_csv(run / "compare" / f"{sid}_{pol}_selection_categories.csv", index_col=0)
    ft = pd.read_csv(next((run / "mzmine_out").glob("*_full_feature_table.csv")), low_memory=False,
                     usecols=["id", "height"]).set_index("id")
    cats = cats.join(ft)
    mgf = OUT_BIG_ROOT / sid / f"{sid}_{pol}_sirius.mgf"

    # MS2 가 있는 봉우리와 그 블록들(문자열 그대로)
    blocks, has_ms2, cur = {}, set(), []
    with open(mgf, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            cur.append(line)
            if line.startswith("END IONS"):
                txt = "".join(cur); cur = []
                fid = next((l.split("=", 1)[1].strip() for l in txt.splitlines() if l.startswith("FEATURE_ID=")), None)
                if fid is None:
                    continue
                blocks.setdefault(int(fid), []).append(txt)
                if "\nMSLEVEL=2" in txt:
                    has_ms2.add(int(fid))
    cats["has_ms2"] = cats.index.isin(has_ms2)
    cats["in_mass_range"] = cats["mz"].between(pr["mz_min"], pr["mz_max"])
    ok = cats[cats.has_ms2 & cats.in_mass_range]

    site = ok[ok.category.str.endswith(("_only", "_dominant"))]
    common = ok[ok.category.isin(["common_similar", "common_present"]) & (ok.qc_rsd <= QC_RSD_MAX)]
    common = common.sort_values("height", ascending=False).head(args.common_top)
    sel = pd.concat([site, common]).sort_values("height", ascending=False)
    sel["target_reason"] = sel["category"]

    out = run / "targets"
    out.mkdir(parents=True, exist_ok=True)
    sel.to_csv(out / f"{sid}_targets.csv", encoding="utf-8-sig", index_label="feature_id")
    with (out / f"{sid}_targets.mgf").open("w", encoding="utf-8", newline="\n") as fh:
        for fid in sel.index:
            for txt in blocks[fid]:
                fh.write(txt)
    info = {"mgf_source": str(mgf), "common_top": args.common_top, "qc_rsd_max_common": QC_RSD_MAX,
            "mz_range": [pr["mz_min"], pr["mz_max"]],
            "available_with_ms2_by_category": ok.category.value_counts().to_dict(),
            "selected_by_category": sel.category.value_counts().to_dict(), "selected_total": int(len(sel))}
    (out / f"{sid}_targets.json").write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(info, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
