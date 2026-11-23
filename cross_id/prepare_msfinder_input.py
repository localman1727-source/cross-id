# -*- coding: utf-8 -*-
"""
MS-FINDER 입력 만들기 — SIRIUS 용으로 골라 둔 봉우리를 MS-FINDER 의 .mat 파일로 옮긴다
=========================================================================================
'왜' MS-FINDER 인가?
  SIRIUS 6.5.4 는 계정 구독이 활성이 아니면 분자식 계산조차 하지 않고 바로 꺼진다(2026-10-08 시험).
  MS-FINDER(RIKEN, MS-DIAL 과 같은 곳)는 계정 없이 이 컴퓨터 안에서만 분자식 추정 → 구조 후보
  순위 매기기를 한다. 구조 후보는 프로그램에 딸린 공개 구조 목록(UNPD·KNApSAcK·COCONUT 등,
  로컬 파일)에서만 고르므로, SIRIUS 처럼 '목록에 없는 구조의 지문'을 예측하지는 못한다.

'왜' 첨가 이온을 [M+H]+ 로 고정하나? (음이온 연구는 [M-H]-, 2026-11-06)
  mzmine 이 붙인 ADDUCT 는 상관 신호로 추정한 것이라 틀릴 수 있다(예: 287.0558 을 [M+NH4]+ 로
  붙였지만 [M+H]+ 로 보면 플라보놀 분자식과 맞는다). 첨가 이온을 바꿔 가며 다 돌리면 계산이 몇 배가
  되므로 1차는 양이온에서 가장 흔한 [M+H]+ 로 통일하고, mzmine 의 추정은 표에 남겨 나중에 대조한다.

'왜' MS1 을 전구 이온 근처만 남기나?
  mgf 의 MS1 블록은 '같이 움직이는 신호' 전체(수십~수백 개)다. MS-FINDER 는 MS1 을 동위원소
  비율 확인에 쓰므로, 전구 이온 −0.5 ~ +3.5 Da 만 남긴다.

입력  runs/<연구>/sirius/<연구>_sirius_input.mgf          (prepare_sirius_input.py)
산출  C:/MSFINDER/work/<연구>/input/<FEATURE_ID>.mat · msfinder_method.ini   (작업 파일 — Dropbox 밖)
      runs/<연구>/msfinder/input_index.csv · msfinder_method.ini (기록용 사본)
'왜' Dropbox 밖인가? MS-FINDER 는 봉우리마다 폴더·파일을 만들고 지운다. Dropbox 가 그 사이 파일을
  붙잡아 지우기 실패가 났다(10-08). mzmine 큰 산출물을 C:/mzmine_out 에 두는 것과 같은 이유.
실행  py scripts/prepare_msfinder_input.py --config config/<연구>.json [--targets] [--spectrum last|consensus] [--limit N]
      --targets   입력을 runs/<연구>/targets/<연구>_targets.mgf(선발 범주 기준, select_annotation_targets.py)로
      --spectrum  last = 봉우리당 마지막 MS2(첫 실행 10-08 의 규칙) / consensus = 런 합본(ms2_consensus.py)
      작업 폴더는 C:/MSFINDER/work/<연구>/<tag>/ (tag = [targets_]<spectrum>; 옛 기본 실행은 tag 없이 <연구>/ 바로 아래)
      ⚠️ 10-08 첫 실행은 last 규칙이었다 — 봉우리마다 런별 MS2 가 최대 12개인데 마지막(대개 QC_03) 하나만 썼다.
"""
import sys
import csv
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console
from ms2_consensus import read_ms2_blocks, consensus

utf8_console()
WORK_ROOT = Path(r"C:\MSFINDER\work")

# MS-FINDER 설정값 — 기본 MSFINDER.INI(v3.73)에서 바꾼 것만 적는다. 나머지는 프로그램 기본값.
METHOD_OVERRIDES = {
    "Ms1Tolerance": "0.01",          # Da. Waters QTof 질량 정확도(수 ppm)보다 넉넉하게
    "MassToleranceType": "Da",
    "Ms2Tolerance": "0.025",         # Da. 조각 이온은 MS1 보다 덜 정확
    "StructureMaximumReportNumber": "10",
    "FormulaMaximumReportNumber": "20",
    "MinesNeverUse": "True",         # MINEs = 반응 규칙으로 '만든' 가상 대사체 → 쓰지 않음
    "MinesOnlyUseForNecessary": "False",
    "MinesAllways": "False",
    "PubChemNeverUse": "True",       # 인터넷 조회 안 함(컴퓨터 안에서만)
    "PubChemOnlyUseForNecessary": "False",
    "PubChemAllways": "False",
    "AllProcess": "True",            # 분자식 + 구조 둘 다
    "FormulaFinder": "True",
    "StructureFinder": "True",
    "TryTopNMolecularFormulaSearch": "5",
}


def read_blocks(path):
    """mgf → [{헤더 dict, peaks [(mz, int)]}]."""
    out, head, peaks, inside = [], {}, [], False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s == "BEGIN IONS":
                head, peaks, inside = {}, [], True
            elif s == "END IONS":
                out.append({"head": head, "peaks": peaks})
                inside = False
            elif inside and s:
                if s[0].isdigit():
                    a = s.split()
                    peaks.append((float(a[0]), float(a[1])))
                elif "=" in s:
                    k, v = s.split("=", 1)
                    head[k] = v
    return out


def write_method(src_ini, dst_ini):
    lines = Path(src_ini).read_text(encoding="utf-8", errors="replace").splitlines()
    seen = set()
    new = []
    for ln in lines:
        k = ln.split("=", 1)[0].strip() if "=" in ln else None
        if k in METHOD_OVERRIDES:
            new.append(f"{k}={METHOD_OVERRIDES[k]}")
            seen.add(k)
        else:
            new.append(ln)
    missing = set(METHOD_OVERRIDES) - seen
    if missing:
        raise SystemExit(f"[중단] 기본 INI 에 없는 키: {sorted(missing)} — MS-FINDER 판이 바뀌었는지 확인")
    Path(dst_ini).write_text("\r\n".join(new) + "\r\n", encoding="utf-8", newline="")


PRECURSOR = {"POS": ("[M+H]+", "Positive"), "NEG": ("[M-H]-", "Negative")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--targets", action="store_true")
    ap.add_argument("--spectrum", choices=["last", "consensus"], default="last")
    ap.add_argument("--limit", type=int, default=None, help="시험용: 앞에서 N 개만")
    ap.add_argument("--msfinder-ini", default=r"C:\MSFINDER\v3.73\MSFINDER v3.73\MSFINDER.INI")
    args = ap.parse_args()

    cfg = load_study(args.config)
    run, sid = cfg["run_dir"], cfg["study_id"]
    # 2026-11-06: 극성별 고정 첨가 이온 — 음이온은 [M-H]- (양이온 고정값을 그대로 쓰면 에러 없이 분자식이 틀린다)
    ptype, imode = PRECURSOR[cfg["polarity"].upper()]
    mgf = run / "targets" / f"{sid}_targets.mgf" if args.targets else run / "sirius" / f"{sid}_sirius_input.mgf"
    tag = ("targets_" if args.targets else "") + args.spectrum
    default_run = (not args.targets) and args.spectrum == "last"          # 10-08 첫 실행과 같은 자리 유지
    out = run / "msfinder" if default_run else run / "msfinder" / tag
    work = WORK_ROOT / sid if default_run else WORK_ROOT / sid / tag
    inp = work / "input"
    out.mkdir(parents=True, exist_ok=True)
    inp.mkdir(parents=True, exist_ok=True)
    for old in inp.glob("*.mat"):
        old.unlink()                              # 지난 입력이 섞이지 않게(우리가 만든 파일만 있는 폴더)

    feats = {}
    for b in read_blocks(mgf):
        fid = b["head"].get("FEATURE_ID")
        lvl = b["head"].get("MSLEVEL")
        feats.setdefault(fid, {})[lvl] = b

    cons = read_ms2_blocks(mgf) if args.spectrum == "consensus" else {}
    order = list(feats)                           # mgf 순서 = 신호 큰 순서
    if args.limit:
        order = order[: args.limit]

    idx = []
    for fid in order:
        f = feats[fid]
        if "2" not in f:
            continue
        h2 = f["2"]["head"]
        pmz = float(h2["PEPMASS"].split()[0])
        rt_min = float(h2["RTINSECONDS"]) / 60.0
        ms1 = [(m, i) for m, i in f.get("1", {"peaks": []})["peaks"] if pmz - 0.5 <= m <= pmz + 3.5]
        ms2 = [(m, i) for m, i in f["2"]["peaks"] if m <= pmz + 0.5]
        cinfo = {}
        if args.spectrum == "consensus" and fid in cons:
            cp, cinfo = consensus(cons[fid])
            ms2 = [(round(m, 5), round(i * 10000, 1)) for m, i in cp]   # 가장 큰 조각 = 10000 (MS-FINDER 는 상대 세기만 본다)
        if not ms2:
            continue
        L = [f"NAME: {fid}", "SCANNUMBER: -1", f"RETENTIONTIME: {rt_min:.4f}",
             f"PRECURSORMZ: {pmz}", f"PRECURSORTYPE: {ptype}", f"IONMODE: {imode}",
             "SPECTRUMTYPE: Centroid", f"INTENSITY: {h2.get('FEATURE_MS1_HEIGHT', '0')}",
             "INCHIKEY: ", "SMILES: ", "FORMULA: ",
             f"COMMENT: mzmine_adduct={h2.get('ADDUCT', '')}",
             "MSTYPE: MS1", f"Num Peaks: {len(ms1)}"]
        L += [f"{m}\t{i}" for m, i in ms1]
        L += ["MSTYPE: MS2", f"Num Peaks: {len(ms2)}"]
        L += [f"{m}\t{i}" for m, i in ms2]
        # newline="" : 윈도우가 \r 을 한 번 더 붙이면(\r\r\n) MS-FINDER 가 빈 줄을 읽다 멈춘다(10-08 실측)
        (inp / f"{fid}.mat").write_text("\r\n".join(L) + "\r\n", encoding="utf-8", newline="")
        idx.append({"feature_id": fid, "precursor_mz": pmz, "rt_min": round(rt_min, 4),
                    "mzmine_adduct": h2.get("ADDUCT", ""), "n_ms1_peaks": len(ms1), "n_ms2_peaks": len(ms2),
                    "collision_energy": h2.get("COLLISION_ENERGY", ""),
                    "spectrum": args.spectrum, "consensus_mode": cinfo.get("mode", ""),
                    "consensus_n_used": cinfo.get("n_used", ""), "n_spectra": cinfo.get("n_spectra", "")})

    with (out / "input_index.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(idx[0]))
        w.writeheader()
        w.writerows(idx)
    write_method(args.msfinder_ini, work / "msfinder_method.ini")
    (out / "msfinder_method.ini").write_bytes((work / "msfinder_method.ini").read_bytes())
    print(f"[생성] {len(idx)} 개 .mat → {inp}")
    print(f"[생성] {work / 'msfinder_method.ini'} (사본: {out})")


if __name__ == "__main__":
    main()
