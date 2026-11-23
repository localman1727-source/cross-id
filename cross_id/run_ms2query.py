# -*- coding: utf-8 -*-
"""
MS2Query 실행 — 공개 스펙트럼 사전에서 '같은 물질' 또는 '가장 닮은 친척 물질'을 찾는다
=========================================================================================
'왜' MS2Query 인가?
  SIRIUS 는 라이선스 승인 전이라 돌지 않고(2026-10-08), MS-FINDER 는 구조 목록 안에서만 고른다.
  MS2Query 는 조각무늬를 학습된 모델(MS2DeepScore·Spec2Vec)로 비교해 공개 사전(GNPS 등에서 정리한
  스펙트럼 묶음)에서 똑같은 물질이 없어도 '가장 닮은 물질'을 돌려준다. 이 컴퓨터 안에서만 돈다.
  ⚠️ 사전에 친척조차 없으면 엉뚱한 이름이 낮은 점수로 붙는다 — 점수(ms2query_model_prediction)와
     전구 이온 질량 차(precursor_mz_difference)를 같이 봐야 한다. 질량 차가 0 근처면 '같은 물질
     후보', 크면 '그만큼 무엇이 붙거나 떨어진 친척'이다.

'왜' 무늬 고르는 방식이 둘인가? (--spectrum)
  sirius_input.mgf 에는 봉우리마다 런별 MS2 가 최대 12개 들어 있다.
  last      = 마지막 MS2 블록 하나(대개 QC_03 런). 첫 MS-FINDER 실행(10-08)과 같은 규칙 — 두 도구를 맞댈 때 쓴다.
  consensus = 여러 런에서 반복된 조각만 남긴 합본(ms2_consensus.py). 잡음이 줄어 이름이 더 잘 붙는지
              last 와 나란히 돌려 결과로 판단한다(10-08 결정).

모델·사전  C:/MS2Query/models_pos  (download_zenodo_files('positive'), Zenodo 13348638)
           C:/MS2Query/models_neg  (download_zenodo_files('negative')) — 설정 polarity 로 고른다(2026-11-06, 음이온 연구 이전 연구)
입력  runs/<연구>/sirius/<연구>_sirius_input.mgf
작업  C:/MS2Query/work/<연구>/<spectrum>/input.mgf · results/   (Dropbox 밖 — 큰 임시 파일)
산출  runs/<연구>/ms2query/<연구>_ms2query_<spectrum>_results.csv · _<spectrum>_input_index.csv
실행  "%LOCALAPPDATA%/venvs/ms2query/Scripts/python.exe" scripts/run_ms2query.py --config config/<연구>.json --spectrum last|consensus [--targets] [--limit N]
      --targets: 입력을 runs/<연구>/targets/<연구>_targets.mgf(선발 범주 기준)로, 산출 이름에 targets_ 가 붙는다
      ('왜' 별도 python? MS2Query 1.5.4 는 torch<2.6·onnxruntime<1.16 을 요구해 Python 3.11 전용 venv 가 필요)
"""
import sys
import csv
import shutil
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console
from ms2_consensus import read_ms2_blocks, consensus

utf8_console()
MODEL_DIRS = {"POS": Path(r"C:\MS2Query\models_pos"), "NEG": Path(r"C:\MS2Query\models_neg")}
# '왜' 극성별인가? 음이온 무늬를 양이온 모델에 넣으면 에러 없이 엉뚱한 이름이 붙는다(2026-11-06)
ION = {"POS": ("positive", "1+"), "NEG": ("negative", "1-")}
WORK_ROOT = Path(r"C:\MS2Query\work")
N_ANALOGS = 3          # 봉우리당 저장할 후보 수(1순위 + 대안 2개)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--spectrum", choices=["last", "consensus"], required=True)
    ap.add_argument("--targets", action="store_true",
                    help="대상을 선발 범주 목록(select_annotation_targets.py)으로 — 산출 이름에 targets_ 가 붙는다")
    ap.add_argument("--limit", type=int, default=None, help="시험용: 앞에서 N 개만")
    ap.add_argument("--mgf", default=None, help="입력 mgf 를 직접 지정(예: 두 산지 차이 봉우리 중 구조 후보가 없던 것). --tag 와 함께")
    ap.add_argument("--tag", default=None, help="--mgf 와 함께: 산출 이름 꼬리표(기존 산출을 덮지 않게)")
    args = ap.parse_args()

    cfg = load_study(args.config)
    run, sid = cfg["run_dir"], cfg["study_id"]
    pol = cfg["polarity"].upper()
    model_dir = MODEL_DIRS[pol]
    # 2026-11-06: 내려받기가 중간에 끊겨 5 MB 짜리 sqlite 만 남은 적이 있다 → 파일 여덟 개가 다 있고 사전이 100 MB 넘는지 본다
    need = ["ms2query_library.sqlite", "ms2ds_embeddings.parquet", "ms2deepscore_model.pt", "spec2vec_model.model",
            "spec2vec_model.model.wv.vectors.npy", "spec2vec_model.model.syn1neg.npy", "ms2query_model.onnx", "s2v_embeddings.parquet"]
    lack = [n for n in need if not (model_dir / n).exists()]
    if lack or (model_dir / need[0]).stat().st_size < 100e6:
        raise SystemExit(f"[중단] {pol} 모델이 없거나 덜 받아졌습니다: {model_dir} (없음 {lack}) — md5 확인 후 다시 받기")
    ionmode, charge = ION[pol]
    mgf = run / "targets" / f"{sid}_targets.mgf" if args.targets else run / "sirius" / f"{sid}_sirius_input.mgf"
    tag = ("targets_" if args.targets else "") + args.spectrum
    if args.mgf:                                   # 2026-10-19 T19: 기본 동작은 그대로, 지정할 때만 바뀐다
        if not args.tag:
            raise SystemExit("[중단] --mgf 를 쓰면 --tag 도 주세요 (기존 산출 덮어쓰기 방지)")
        mgf, tag = Path(args.mgf), f"{args.tag}_{args.spectrum}"
    work = WORK_ROOT / sid / tag
    work.mkdir(parents=True, exist_ok=True)
    res_dir = work / "results"
    if res_dir.exists():
        shutil.rmtree(res_dir)                       # 지난 결과 파일 이름 뒤에 번호가 붙어 섞이는 것 방지
    out = run / "ms2query"
    out.mkdir(parents=True, exist_ok=True)

    feats = read_ms2_blocks(mgf)
    order = list(feats)[: args.limit] if args.limit else list(feats)
    L, idx = [], []
    for fid in order:
        sp = feats[fid]
        if args.spectrum == "last":
            peaks, info = sp[-1]["peaks"], {"n_spectra": len(sp), "n_used": 1, "median_cosine": None, "mode": "last",
                                                "used_files": sp[-1]["file"]}
            src = sp[-1]["file"]
        else:
            peaks, info = consensus(sp)
            src = "consensus"
        f = sp[-1]
        # 머리줄은 matchms 가 읽는 표준 이름만: feature_id 는 결과표에 다시 붙이려고 남긴다
        L += ["BEGIN IONS", f"FEATURE_ID={fid}", f"PEPMASS={f['pmz']}", f"CHARGE={charge}",
              f"IONMODE={ionmode}", f"RTINSECONDS={f['rt_sec']}"] + [f"{m} {i}" for m, i in peaks] + ["END IONS", ""]
        idx.append({"feature_id": fid, "precursor_mz": f["pmz"], "rt_min": round(f["rt_sec"] / 60, 4),
                    "ms2_source": src, "mzmine_adduct": f["adduct"], "n_ms2_peaks": len(peaks), **info})
    (work / "input.mgf").write_text("\n".join(L), encoding="utf-8")
    with (out / f"{sid}_ms2query_{tag}_input_index.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(idx[0]))
        w.writeheader()
        w.writerows(idx)
    print(f"[입력] 봉우리 {len(idx)} 개 → {work / 'input.mgf'}")

    from ms2query.ms2library import create_library_object_from_one_dir
    from ms2query.run_ms2query import run_ms2query_single_file
    from ms2query.utils import SettingsRunMS2Query
    lib = create_library_object_from_one_dir(str(model_dir))
    settings = SettingsRunMS2Query(nr_of_top_analogs_to_save=N_ANALOGS,
                                   additional_metadata_columns=("feature_id", "retention_time"),
                                   additional_ms2query_score_columns=("s2v_score", "ms2ds_score"))
    run_ms2query_single_file("input.mgf", str(work), str(res_dir), lib, settings)
    res = next(res_dir.glob("input*.csv"))
    shutil.copyfile(res, out / f"{sid}_ms2query_{tag}_results.csv")
    print(f"[생성] {out / f'{sid}_ms2query_{tag}_results.csv'}")


if __name__ == "__main__":
    main()
