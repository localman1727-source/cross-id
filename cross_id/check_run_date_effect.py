# -*- coding: utf-8 -*-
"""
S2 측정일 효과 점검 — 같은 시료를 다른 날 잰 차이가 산지 차이만큼 큰가?
========================================================================
'왜' 필요한가? (2026-11-06 결정 Q1, ANALYSIS_DESIGN_2026-11-06.md 4절)
  이전 연구 은 시료마다 1반복을 8/14 에, 2·3반복을 8/25 에 연달아 주입했다. 세 런을 모두 반복으로
  쓰기로 했으므로, 측정일이 바뀐 것만으로 생기는 차이가 산지 사이 차이와 같은 크기라면
  "산지 차이"의 일부가 날짜 차이일 수 있다 → 그 경우 결정을 다시 연다.

비교하는 세 가지 차이 (모두 log2 면적, 시료 모든 런에서 잡힌 봉우리만)
  date_shift      같은 시료: 혼자 잰 날 런 − 같은 날 연달아 잰 두 런 평균
  same_day        같은 시료: 같은 날 연달아 잰 두 런끼리 차이 (가장 작은 흔들림 기준)
  site_sameday    다른 시료끼리: 같은 날 런 평균의 차이 (날짜가 섞이지 않은 산지 차이)
  추가로 날짜 차이가 네 시료 모두 같은 방향으로 움직이는지(기기 쪽 체계적 이동인지) 센다.

QC 는 8/25 에만 있어 날짜 비교에 쓰지 않는다. compare_samples.py 결과는 건드리지 않는다.
날짜별 검출 수(범주 표에 붙일 칸)는 <sid>_detect_by_date.csv 로 따로 남긴다 — 봉우리 id 로 이어 붙인다.

쓰는 법:  py scripts\\check_run_date_effect.py --config config\\<연구>.json
산출물:  runs/<연구>/compare/<sid>_run_date_check.{json,md} · _date_effect_features.csv · _detect_by_date.csv
"""
import sys
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import utf8_console, add_config_arg, load_study

FOLD2_LOG2 = 1.0     # 두 배 차이 — compare_samples.FOLD_MIN_LOG2 와 같은 선


def main():
    utf8_console()
    ap = add_config_arg(argparse.ArgumentParser())
    ap.parse_args()
    cfg = load_study()
    run_dir = cfg["run_dir"]
    sid = f"{cfg['study_id']}_{cfg['polarity']}"          # compare_samples.py 와 같은 이름 규칙
    out = run_dir / "compare"
    out.mkdir(exist_ok=True)

    man = pd.read_csv(run_dir / "raw_manifest.csv", encoding="utf-8-sig")
    # '측정일' = 묶음(배치) 날짜 = 런 이름 앞 8자리. 획득 시각으로 묶으면 자정을 넘긴 연속 주입이
    # 다른 날로 갈라진다(2026-11-06 이전 연구 한 시료 런이 8/26 00시대 획득 — 이름은 20260825 묶음).
    man["date"] = man["run_name"].str[:8]
    if not man["date"].str.fullmatch(r"\d{8}").all():
        raise SystemExit("[중단] 런 이름 앞 8자리가 날짜(YYYYMMDD)가 아닌 런이 있습니다 — 묶음 날짜 규칙을 확인하세요")
    man["order"] = pd.to_datetime(man["acquired_date"] + " " + man["acquired_time"]).rank().astype(int)
    ft = pd.read_csv(run_dir / "mzmine_out" / f"{sid}_full_feature_table.csv", low_memory=False)
    col = {r.run_name: f"datafile:{r.run_name}.raw:area" for r in man.itertuples()}
    miss = [c for c in col.values() if c not in ft.columns]
    if miss:
        raise SystemExit(f"[중단] 봉우리 표에 없는 런: {miss[:3]}")
    X = ft[list(col.values())].copy()
    X.columns = list(col.keys())
    X.index = ft["id"]

    sm = man[man.is_qc == "N"]
    samples = list(sm.sample_code.unique())
    # 시료마다 날짜별 런 묶기: '혼자 잰 날'(런 1개)과 '연달아 잰 날'(런 2개 이상)
    design = {}
    for s in samples:
        g = sm[sm.sample_code == s].sort_values("order").groupby("date")["run_name"].apply(list).to_dict()
        single = [d for d, r in g.items() if len(r) == 1]
        multi = [d for d, r in g.items() if len(r) >= 2]
        if len(single) != 1 or len(multi) != 1:
            raise SystemExit(f"[중단] {s}: 날짜 구성이 '하루 1런 + 하루 2런 이상'이 아닙니다 {g}")
        design[s] = {"single_date": single[0], "single_run": g[single[0]][0],
                     "multi_date": multi[0], "multi_runs": g[multi[0]]}

    # 날짜별 검출 수 — 전체 봉우리 표 기준(범주 표와 같은 '검출' = 빈칸 아님 & 0 초과)
    det = (X.notna() & (X > 0))
    dbd = pd.DataFrame(index=X.index)
    for s, d in design.items():
        dbd[f"{s}_n_{d['single_date']}"] = det[d["single_run"]].astype(int)
        dbd[f"{s}_n_{d['multi_date']}"] = det[d["multi_runs"]].sum(axis=1)
    dbd.to_csv(out / f"{sid}_detect_by_date.csv", encoding="utf-8-sig")

    # 이후는 모든 시료 런에서 잡힌 봉우리만 (빈칸이 차이를 부풀리지 않게)
    runs = sm.run_name.tolist()
    full = det[runs].all(axis=1)
    L = np.log2(X.loc[full, runs])
    f = pd.DataFrame(index=L.index)
    f["mz"] = ft.set_index("id").loc[L.index, "mz"]
    f["rt"] = ft.set_index("id").loc[L.index, "rt"]
    for s, d in design.items():
        m = L[d["multi_runs"]]
        f[f"{s}_date_shift"] = L[d["single_run"]] - m.mean(axis=1)
        f[f"{s}_same_day"] = m.iloc[:, 0] - m.iloc[:, 1] if m.shape[1] == 2 else m.max(axis=1) - m.min(axis=1)
        f[f"{s}_multi_mean"] = m.mean(axis=1)
    site_cols = []
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            c = f"site_sameday_{a}_vs_{b}"
            f[c] = f[f"{a}_multi_mean"] - f[f"{b}_multi_mean"]
            site_cols.append(c)
    ds = f[[f"{s}_date_shift" for s in samples]]
    f["date_shift_median_abs"] = ds.abs().median(axis=1)
    f["date_shift_same_sign_all"] = (np.sign(ds).nunique(axis=1) == 1) & (ds.abs() > 0).all(axis=1)
    f["date_shift_mean"] = ds.mean(axis=1)
    f["site_sameday_median_abs"] = f[site_cols].abs().median(axis=1)
    f.to_csv(out / f"{sid}_date_effect_features.csv", encoding="utf-8-sig")

    def summ(v):
        v = v.abs().stack() if isinstance(v, pd.DataFrame) else v.abs()
        return {"median_abs_log2": float(v.median()), "p90_abs_log2": float(v.quantile(0.9)),
                "frac_ge_2fold": float((v >= FOLD2_LOG2).mean()), "n_values": int(v.size)}

    sd = f[[f"{s}_same_day" for s in samples]]
    res = {
        "study": sid, "n_features_all": int(len(X)), "n_features_all_sample_runs": int(full.sum()),
        "design": design, "threshold_fold_log2": FOLD2_LOG2,
        "date_shift": summ(ds), "same_day": summ(sd), "site_sameday": summ(f[site_cols]),
        "date_shift_per_sample": {s: summ(f[f"{s}_date_shift"]) for s in samples},
        "same_day_per_sample": {s: summ(f[f"{s}_same_day"]) for s in samples},
        "site_sameday_per_pair": {c: summ(f[c]) for c in site_cols},
        "date_shift_systematic": {
            "frac_same_sign_all_samples": float(f["date_shift_same_sign_all"].mean()),
            "frac_same_sign_and_mean_ge_2fold": float((f["date_shift_same_sign_all"] & (f["date_shift_mean"].abs() >= FOLD2_LOG2)).mean()),
            "median_of_feature_mean_shift_log2": float(f["date_shift_mean"].median()),
        },
    }
    # 런 사이 순위 상관: 같은 시료 다른 날 / 같은 시료 같은 날 / 다른 시료 같은 날
    C = L.corr(method="spearman")
    same_diff, same_same, other_same = [], [], []
    for s, d in design.items():
        for r in d["multi_runs"]:
            same_diff.append(C.loc[d["single_run"], r])
        mr = d["multi_runs"]
        same_same += [C.loc[mr[i], mr[j]] for i in range(len(mr)) for j in range(i + 1, len(mr))]
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            same_same_ab = [C.loc[x, y] for x in design[a]["multi_runs"] for y in design[b]["multi_runs"]]
            other_same += same_same_ab
    res["run_rank_correlation"] = {
        "same_sample_diff_day": {"median": float(np.median(same_diff)), "min": float(np.min(same_diff)), "n": len(same_diff)},
        "same_sample_same_day": {"median": float(np.median(same_same)), "min": float(np.min(same_same)), "n": len(same_same)},
        "diff_sample_same_day": {"median": float(np.median(other_same)), "max": float(np.max(other_same)), "n": len(other_same)},
    }
    r = res["date_shift"]["median_abs_log2"] / res["site_sameday"]["median_abs_log2"]
    res["ratio_date_over_site_median_abs"] = float(r)
    res["ratio_date_over_sameday_median_abs"] = float(res["date_shift"]["median_abs_log2"] / res["same_day"]["median_abs_log2"])
    # 판정 문구: 설계 4절 "측정일 차이가 산지 차이와 같은 크기면 재논의" — 같은 크기 = 비율 0.5 이상으로 둔다(기록용 선)
    res["reopen_rule"] = "date_shift median|log2| >= 0.5 x site_sameday median|log2| → Q1 재논의"
    res["reopen_flag"] = bool(r >= 0.5)
    (out / f"{sid}_run_date_check.json").write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [f"# 측정일 효과 점검 — {sid}", "",
          f"봉우리 전체 {res['n_features_all']} · 모든 시료 런 검출 {res['n_features_all_sample_runs']}", "",
          "| 차이 | median abs log2 | p90 abs log2 | 2배 이상 비율 | n |", "|---|---|---|---|---|"]
    for k in ("same_day", "date_shift", "site_sameday"):
        v = res[k]
        md.append(f"| {k} | {v['median_abs_log2']} | {v['p90_abs_log2']} | {v['frac_ge_2fold']} | {v['n_values']} |")
    md += ["", f"date/site 비율 {r} · date/same_day 비율 {res['ratio_date_over_sameday_median_abs']} · 재논의 표시 {res['reopen_flag']}",
           "", "체계적 이동: " + json.dumps(res["date_shift_systematic"], ensure_ascii=False),
           "", "런 순위 상관: " + json.dumps(res["run_rank_correlation"], ensure_ascii=False), "",
           "시료별 date_shift: " + json.dumps(res["date_shift_per_sample"], ensure_ascii=False), ""]
    (out / f"{sid}_run_date_check.md").write_text("\n".join(md), encoding="utf-8")
    print(f"[생성] {out / (sid + '_run_date_check.json')} (값은 파일에)")


if __name__ == "__main__":
    main()
