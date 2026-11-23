# -*- coding: utf-8 -*-
"""
⑦ 시료 비교 — 같은 시료 반복의 흩어짐 vs 시료 사이 차이
=========================================================
mzmine 전체 봉우리 표(full_feature_table)와 이름 후보 표(annotations)를 읽어
  1) 반복확인용(QC) 시료로 안정한 봉우리만 남기고
  2) 봉우리마다 시료 평균과 반복 흩어짐을 비교해 "반복보다 큰 차이"를 골라내고
  3) 주입 순서에 따른 흔들림이 그 차이와 같은 크기인지 점검하고
  4) 이름 후보 중 성립하기 어려운 이온 형태를 걸러 후보 표를 만든다.
  5) 선발 범주를 매긴다(2026-10-08 추가) — 공통 / 산지 특이 / 산지 우세.

'왜' 5) 는 1) 의 '모든 런에서 잡힘' 조건을 쓰지 않나?
  1)~3) 은 모든 런에서 잡힌 봉우리만 다루므로, 한 산지에만 있는 물질은 처음부터 빠진다.
  산지 특이 물질을 고르려면 '다른 산지에서 안 잡힘' 자체가 정보이므로 전체 봉우리 표에서 따로 본다.
    common_similar  세 산지 반복 9런 모두 검출 + 산지 평균 최대/최소 < 2배(FOLD_MIN_LOG2)
                    → 종의 공통 성분 후보(산지와 무관한 기준선)
    common_present  세 산지 9런 모두 검출(양 차이는 따지지 않음; common_similar 를 포함하지 않는 나머지)
    <산지>_only     그 산지 3런 모두 검출 + 다른 두 산지 6런 모두 미검출
    <산지>_dominant 그 산지 3런 모두 검출 + 다른 두 산지 각각보다 10배(SITE_DOMINANT_LOG2) 이상 많고
                    반복 범위가 겹치지 않음(그 산지 최솟값 > 다른 산지 최댓값)
  ⚠️ '미검출' = mzmine 의 빈칸 메우기(gap filling) 뒤에도 빈 칸. 아주 적게 있는 것을 못 잡았을 수 있다.
  ⚠️ 산지당 채집 시료가 하나라 '산지 특이'는 그 한 그루·그 채집 시점의 특이일 수 있다(연도·보관 차이 포함).
  QC 는 범주 조건에 넣지 않고 칸(qc_detected, qc_rsd)으로만 남긴다 — QC 가 세 산지를 섞은 것인지
  확인 전이라, 한 산지 물질이 QC 에서 약하게 잡히는 것을 탈락 사유로 삼지 않는다.

'왜' 검정 수치(p)를 내지 않나?
  반복의 종류는 연구마다 다르므로 설정의 replicate_type 에 적는다. 첫 연구는
  같은 채집 시료를 따로 세 번 추출해 한 번씩 주입한 추출 반복이다(2026-10-08 확인 —
  처음엔 한 추출물을 세 번 주입한 기술 반복으로 잘못 적었다). 추출 반복은 추출·주입 흩어짐을
  함께 담지만, 산지마다 채집 시료가 하나(n=1)인 것은 그대로다. 반복 셋으로 낸 p 는
  "같은 시료를 세 번 추출한 일치도"이지 산지 차이의 확실성이 아니다(pseudo-replication).
  그래서 방향·크기·반복 흩어짐과의 비교만 남긴다(PROCESS.md ⑨).

'왜' 이온 형태를 거르나?
  후보 대조는 질량 하나만 맞추는 방식이라, 작은 분자에 전하 두 개([M+2H]2+)나
  전자가 빠진 형태([M-e]+)를 붙여도 질량이 맞으면 이름이 붙는다. 양이온 전기분무에서
  이런 형태가 작은 천연물의 주 이온일 가능성은 낮다. 지우지 않고 등급을 낮춘다.

쓰는 법:
    py scripts\\compare_samples.py --config config\\<연구>.json
산출물: runs/<연구>/compare/
"""
import sys
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import utf8_console, add_config_arg, load_study

# ── 판정 기준 (결과 JSON 에 함께 기록된다) ────────────────────────────
QC_RSD_MAX = 0.30          # QC 세 번의 변동계수 상한 — 대사체 분야 관행값
FOLD_MIN_LOG2 = 1.0        # 두 배 이상 차이
SITE_DOMINANT_LOG2 = float(np.log2(10))   # 산지 우세 = 다른 두 산지 각각보다 10배 이상
PPM_MAX = 5.0              # 이름 후보 질량 오차 상한
# 극성별 이온 형태 (2026-11-06: 음이온 연구 이전 연구 추가). ADDUCT_OK 는 기록용이고 등급은 DIMER·UNLIKELY 로만 매긴다.
#  '왜' 극성별인가? 양이온 집합만 두면 음이온의 이량체·두 전하 형태가 걸리지 않고 B_mass_only 로 한 칸 높게 남는다.
ADDUCTS = {
    "POS": {"ok": {"[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M-H2O+H]+"},
            "dimer": {"[2M+H]+", "[2M+Na]+"},              # 진한 봉우리에서만 흔함 → 낮은 등급
            "unlikely": {"[M+2H]2+", "[M-e]+"}},           # 작은 분자에서 주 이온일 가능성 낮음
    "NEG": {"ok": {"[M-H]-", "[M+HCOO]-", "[M+FA-H]-", "[M+Cl]-", "[M-H2O-H]-"},
            "dimer": {"[2M-H]-", "[2M+HCOO]-", "[2M+FA-H]-"},
            "unlikely": {"[M-2H]2-", "[M+e]-"}},
}
ADDUCT_OK, ADDUCT_DIMER, ADDUCT_UNLIKELY = (ADDUCTS["POS"][k] for k in ("ok", "dimer", "unlikely"))


def site_categories(X, reps, qc_runs, samples, ft, confirm=None):
    """선발 범주(공통·산지 특이·산지 우세). 전체 봉우리 표에서 본다(빈칸 = 미검출).
    confirm: {시료: [확인용 런]} — 범주 조건에는 넣지 않고 칸으로만 남긴다(2026-11-06 이전 연구 Q1 재결정)."""
    det = {s: X[reps[s]].notna() & (X[reps[s]] > 0) for s in samples}
    ndet = pd.DataFrame({s: det[s].sum(axis=1) for s in samples})
    nrep = pd.Series({s: len(reps[s]) for s in samples})
    L = np.log2(X.where(X > 0))
    mean = pd.DataFrame({s: L[reps[s]].mean(axis=1) for s in samples})
    lo = pd.DataFrame({s: L[reps[s]].min(axis=1) for s in samples})
    hi = pd.DataFrame({s: L[reps[s]].max(axis=1) for s in samples})
    full = (ndet == nrep).all(axis=1)
    spread = mean.max(axis=1) - mean.min(axis=1)
    cat = pd.Series("other", index=X.index)
    cat[full] = "common_present"
    cat[full & (spread < FOLD_MIN_LOG2)] = "common_similar"
    for s in samples:
        others = [o for o in samples if o != s]
        only = (ndet[s] == nrep[s]) & (ndet[others] == 0).all(axis=1)
        dom = (ndet[s] == nrep[s]) & ~only
        for o in others:
            # 미검출 산지는 비교 대상에서 '없음'으로 본다(어떤 차이든 성립)
            gap_ok = mean[o].isna() | (mean[s] - mean[o] >= SITE_DOMINANT_LOG2)
            sep_ok = hi[o].isna() | (lo[s] > hi[o])
            dom &= gap_ok & sep_ok
        cat[only] = f"{s}_only"
        cat[dom] = f"{s}_dominant"
    qc = X[qc_runs]
    res = pd.DataFrame({"mz": ft["mz"].values, "rt": ft["rt"].values}, index=X.index)
    res["category"] = cat
    for s in samples:
        res[f"{s}_n_detected"] = ndet[s]
        res[f"{s}_log2mean"] = mean[s]
    res["max_over_min_site_log2"] = spread
    res["qc_detected"] = qc.notna().sum(axis=1)
    res["qc_rsd"] = qc.std(axis=1, ddof=1) / qc.mean(axis=1)
    if confirm:
        # 확인용 런(다른 측정일): 그 시료 확인 런에서 잡혔는지, 산지 범주라면 그 산지 확인 런 모두에서 잡혔는지
        cdet = {s: (X[confirm[s]].notna() & (X[confirm[s]] > 0)).sum(axis=1) for s in samples}
        for s in samples:
            res[f"{s}_n_confirm"] = cdet[s]
        site_of = cat.str.replace(r"_(only|dominant)$", "", regex=True)
        res["site_confirmed"] = pd.Series(pd.NA, index=X.index, dtype="object")
        for s in samples:
            m = site_of == s
            res.loc[m, "site_confirmed"] = (cdet[s][m] == len(confirm[s]))
        others_c = {s: sum(cdet[o] for o in samples if o != s) for s in samples}
        res["others_confirm_detected"] = pd.Series(pd.NA, index=X.index, dtype="object")
        for s in samples:
            m = site_of == s
            res.loc[m, "others_confirm_detected"] = others_c[s][m]
    return res


def main():
    utf8_console()
    ap = add_config_arg(argparse.ArgumentParser())
    ap.parse_args()
    cfg = load_study()
    run_dir = cfg["run_dir"]
    sid = f"{cfg['study_id']}_{cfg['polarity']}"
    global ADDUCT_OK, ADDUCT_DIMER, ADDUCT_UNLIKELY
    ADDUCT_OK, ADDUCT_DIMER, ADDUCT_UNLIKELY = (ADDUCTS[cfg["polarity"].upper()][k] for k in ("ok", "dimer", "unlikely"))
    mz_dir = run_dir / "mzmine_out"
    out = run_dir / "compare"
    out.mkdir(exist_ok=True)

    man = pd.read_csv(run_dir / "raw_manifest.csv", encoding="utf-8-sig")
    man["order"] = pd.to_datetime(man["acquired_date"] + " " + man["acquired_time"]).rank().astype(int)
    ft = pd.read_csv(mz_dir / f"{sid}_full_feature_table.csv", low_memory=False)

    col = {r.run_name: f"datafile:{r.run_name}.raw:area" for r in man.itertuples()}
    missing = [c for c in col.values() if c not in ft.columns]
    if missing:
        raise SystemExit(f"[중단] 봉우리 표에 없는 런이 있습니다 (시험 실행 결과일 수 있음): {missing[:3]}")
    X = ft[list(col.values())].copy()
    X.columns = list(col.keys())
    X.index = ft["id"]

    qc_runs = man.loc[man.is_qc == "Y", "run_name"].tolist()
    samples = [s for s in man.loc[man.is_qc == "N", "sample_code"].unique()]
    reps = {s: man.loc[man.sample_code == s].sort_values("order")["run_name"].tolist() for s in samples}
    # 2026-11-06 이전 연구 Q1 재결정: 설정 replicate_batches(런 이름 앞 8자리 = 묶음 날짜)가 있으면 그 묶음 런만 반복으로,
    # 나머지 시료 런은 확인용 칸으로. 없으면 예전과 같이 시료 런 전부가 반복(첫 연구 등 그대로).
    rb = cfg.get("replicate_batches")
    confirm = None
    if rb:
        confirm = {s: [r for r in reps[s] if r[:8] not in rb] for s in samples}
        reps = {s: [r for r in reps[s] if r[:8] in rb] for s in samples}
        if any(len(v) < 2 for v in reps.values()):
            raise SystemExit(f"[중단] replicate_batches {rb} 로 반복이 2 런 미만인 시료가 있습니다: {reps}")
        print(f"[반복] 묶음 {rb} 만 반복으로 사용, 나머지는 확인용 칸: {sum(map(len, confirm.values()))} 런")

    cats = site_categories(X, reps, qc_runs, samples, ft, confirm)
    cats.to_csv(out / f"{sid}_selection_categories.csv", encoding="utf-8-sig")

    # 1) 모든 런에서 잡힌 봉우리 + QC 안정성
    n_all = len(X)
    use = [r for s in samples for r in reps[s]] + qc_runs
    full = X[use].notna().all(axis=1) & (X[use] > 0).all(axis=1)
    Xf = X.loc[full, use]
    qc = Xf[qc_runs]
    qc_rsd = qc.std(axis=1, ddof=1) / qc.mean(axis=1)
    keep = qc_rsd <= QC_RSD_MAX
    Xk = Xf[keep]
    L = np.log2(Xk)

    # 2) 시료 평균·반복 흩어짐·쌍별 차이
    rows = pd.DataFrame(index=L.index)
    rows["mz"] = ft.set_index("id").loc[L.index, "mz"]
    rows["rt"] = ft.set_index("id").loc[L.index, "rt"]
    rows["qc_rsd"] = qc_rsd[keep]
    for s in samples:
        rows[f"{s}_log2mean"] = L[reps[s]].mean(axis=1)
        rows[f"{s}_log2sd"] = L[reps[s]].std(axis=1, ddof=1)
        rows[f"{s}_min"] = L[reps[s]].min(axis=1)
        rows[f"{s}_max"] = L[reps[s]].max(axis=1)
    pair_summary = {}
    sep_any = pd.Series(False, index=L.index)
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            d = rows[f"{a}_log2mean"] - rows[f"{b}_log2mean"]
            no_overlap = (rows[f"{a}_min"] > rows[f"{b}_max"]) | (rows[f"{b}_min"] > rows[f"{a}_max"])
            sep = (d.abs() >= FOLD_MIN_LOG2) & no_overlap
            rows[f"log2FC_{a}_vs_{b}"] = d
            rows[f"sep_{a}_vs_{b}"] = sep
            sep_any |= sep
            pair_summary[f"{a}_vs_{b}"] = {
                "n_separated": int(sep.sum()),
                "n_higher_in_first": int((sep & (d > 0)).sum()),
                "n_higher_in_second": int((sep & (d < 0)).sum()),
                "median_abs_log2FC_all": float(d.abs().median()),
            }
    rows["separated_any_pair"] = sep_any
    within_sd = pd.concat([rows[f"{s}_log2sd"] for s in samples], axis=1).median(axis=1)
    between_sd = rows[[f"{s}_log2mean" for s in samples]].std(axis=1, ddof=1)
    rows["between_over_within"] = between_sd / within_sd.replace(0, np.nan)

    # 3) 주입 순서 흔들림 점검 — 같은 시료 첫 주입 → 마지막 주입의 차이
    drift = {}
    for s in samples + ["QC"]:
        r = reps[s] if s != "QC" else man.loc[man.is_qc == "Y"].sort_values("order")["run_name"].tolist()
        dd = L[r[-1]] - L[r[0]]
        drift[s] = {"median_log2_last_vs_first": float(dd.median()),
                    "frac_features_abs_ge_1": float((dd.abs() >= 1).mean())}

    # 런 사이 유사도 (전체 봉우리 순위 상관)
    corr = L.corr(method="spearman")
    order_names = [r for r in man.sort_values("order")["run_name"].tolist() if r in L.columns]
    corr.loc[order_names, order_names].to_csv(out / f"{sid}_run_correlation.csv", encoding="utf-8-sig")

    # 주성분 (그림용 좌표만, 평균 중심 + 파레토 크기맞춤)
    Z = L.T
    Z = (Z - Z.mean()) / np.sqrt(Z.std(ddof=1))
    U, S, Vt = np.linalg.svd(Z.values, full_matrices=False)
    ev = (S ** 2) / (S ** 2).sum()
    pcs = pd.DataFrame(U[:, :3] * S[:3], index=Z.index, columns=["PC1", "PC2", "PC3"])
    pcs = pcs.join(man.set_index("run_name")[["sample_code", "is_qc", "order"]])
    pcs.to_csv(out / f"{sid}_pca_scores.csv", encoding="utf-8-sig")

    rows.to_csv(out / f"{sid}_feature_comparison.csv", encoding="utf-8-sig")

    # 4) 이름 후보 정리
    ann = pd.read_csv(mz_dir / f"{sid}_annotations.csv")
    ann = ann[ann["method"].astype(str).str.contains("compound_db", na=False)].copy()
    ftx = ft.set_index("id")
    ann["feature_mz"] = ann["id"].map(ftx["mz"])
    ann["feature_rt"] = ann["id"].map(ftx["rt"])
    ann["ppm"] = (ann["feature_mz"] - ann["precursor_mz"]) / ann["precursor_mz"] * 1e6
    iin = ftx.get("ion_identities:ion_identities")
    ann["iin_adduct"] = ann["id"].map(iin) if iin is not None else np.nan

    def grade(r):
        if abs(r.ppm) > PPM_MAX:
            return "C_mass_off"
        if r.adduct in ADDUCT_UNLIKELY:
            return "C_unlikely_ion"
        if r.adduct in ADDUCT_DIMER:
            return "B_dimer"
        if isinstance(r.iin_adduct, str) and r.adduct in r.iin_adduct:
            return "A_ion_network_agrees"
        return "B_mass_only"
    ann["grade"] = ann.apply(grade, axis=1)
    ann["in_stable_set"] = ann["id"].isin(rows.index)
    ann["separated_any_pair"] = ann["id"].map(rows["separated_any_pair"]).fillna(False)
    for s in samples:
        ann[f"{s}_log2mean"] = ann["id"].map(rows[f"{s}_log2mean"])
    for k in [c for c in rows.columns if c.startswith("log2FC_")]:
        ann[k] = ann["id"].map(rows[k])
    ann.sort_values(["grade", "id"]).to_csv(out / f"{sid}_annotations_graded.csv",
                                            index=False, encoding="utf-8-sig")

    res = {
        "study": sid,
        "criteria": {"qc_rsd_max": QC_RSD_MAX, "fold_min_log2": FOLD_MIN_LOG2,
                     "separation_rule": "|log2FC|>=fold_min AND replicate ranges do not overlap",
                     "ppm_max": PPM_MAX, "adduct_ok": sorted(ADDUCT_OK),
                     "adduct_dimer": sorted(ADDUCT_DIMER), "adduct_unlikely": sorted(ADDUCT_UNLIKELY),
                     "replicate_type": cfg.get("replicate_type", "unspecified"),
                     "no_p_values_reason": "one collected sample per site; replicates are within-sample (pseudo-replication for site inference)"},
        "selection_categories": {
            "replicate_batches": cfg.get("replicate_batches"),
            "confirm_runs": confirm,
            "rule": {"common_similar": "all replicate runs of all sites detected AND max/min site mean < 2-fold",
                     "common_present": "all replicate runs of all sites detected, not common_similar",
                     "<site>_only": "all replicate runs of site detected AND 0 replicate runs of other sites",
                     "<site>_dominant": "all replicate runs of site detected AND >=10-fold over each other site mean "
                                        "(missing = absent) AND site min > other sites' max",
                     "site_dominant_log2": SITE_DOMINANT_LOG2},
            "counts": cats["category"].value_counts().to_dict(),
            "site_confirmed_counts": (cats.loc[cats.category.str.endswith(("_only", "_dominant")), "site_confirmed"].astype(str).value_counts().to_dict() if "site_confirmed" in cats else None)},
        "n_features_total": int(n_all),
        "n_detected_all_runs": int(full.sum()),
        "n_qc_stable": int(keep.sum()),
        "qc_rsd_median_detected_all": float(qc_rsd.median()),
        "replicate_log2sd_median": {s: float(rows[f"{s}_log2sd"].median()) for s in samples},
        "between_over_within_median": float(rows["between_over_within"].median()),
        "pairs": pair_summary,
        "n_separated_any_pair": int(sep_any.sum()),
        "drift_last_vs_first": drift,
        "pca_explained": [float(x) for x in ev[:3]],
        "run_correlation_min_within_sample": {
            s: float(corr.loc[reps[s], reps[s]].where(~np.eye(len(reps[s]), dtype=bool)).min().min())
            for s in samples},
        "run_correlation_max_between_samples": float(max(
            corr.loc[reps[a], reps[b]].values.max()
            for i, a in enumerate(samples) for b in samples[i + 1:])),
        "annotations": {
            "rows_total": int(len(ann)),
            "features_annotated": int(ann["id"].nunique()),
            "by_grade_rows": ann["grade"].value_counts().to_dict(),
            "features_grade_A_or_B_in_stable_set": int(ann.loc[ann.grade.str[0].isin(["A", "B"]) & ann.in_stable_set, "id"].nunique()),
            "features_grade_A_or_B_separated": int(ann.loc[ann.grade.str[0].isin(["A", "B"]) & ann.separated_any_pair, "id"].nunique()),
        },
    }
    (out / f"{sid}_compare_results.json").write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
