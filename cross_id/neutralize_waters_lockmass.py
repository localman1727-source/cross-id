# -*- coding: utf-8 -*-
"""
Waters 원시자료에 '보정 배율 1' lmgt.inf 를 써 넣어 mzmine 의 lock mass 재보정을 막는다
========================================================================
'왜' 필요한가? (2026-11-09 이전 연구)
  측정 단계에서 이미 보정된 원시자료에 mzmine 이 lock mass 보정을 다시 걸어 m/z 가 −20~−80 ppm 밀렸다.
  mzmine 4.10.6 배치의 'Apply lockmass on import (Waters)' 를 false 로 해도 로그에 "Ignoring" 이라 찍힌 뒤
  Waters 읽기 부품이 보정을 계산해 lmgt.inf 를 새로 쓰고 그대로 적용했다(시험: methyl gallate −54.8 ppm).
  lmgt.inf 가 이미 있으면 mzmine 은 "already lock mass corrected" 로 보고 건너뛰므로,
  배율을 모두 1.0 으로 둔 lmgt.inf 를 미리 넣어 두면 보정이 사실상 꺼진다(시험: +2.2 · −2.5 ppm).

동작: raw_manifest 의 모든 .raw 에 대해
  - lmgt.inf 가 없으면 배율 1.0 파일을 쓴다(시간 0~60 분, 1 분 간격).
  - 이미 있고 배율이 전부 1.0 이면 그대로 둔다.
  - 배율 1.0 이 아닌 파일이 있으면 archive_lmgt_written_by_mzmine_<오늘>/ 에 백업한 뒤 1.0 으로 바꾼다.
  --check 는 바꾸지 않고 상태만 센다.
쓰는 법: py scripts/neutralize_waters_lockmass.py --config config/<연구>.json [--check]
"""
import sys
import shutil
import argparse
from datetime import date
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import utf8_console, load_study

PROJECT = Path(__file__).resolve().parents[1]


def is_identity(p: Path) -> bool:
    lines = p.read_text(encoding="ascii", errors="replace").splitlines()[2:]
    vals = [float(l.split(",")[1]) for l in lines if "," in l]
    return bool(vals) and all(abs(v - 1.0) < 1e-9 for v in vals)


def identity_text(lockmass: float) -> str:
    rows = ["1", f"{lockmass:.9f},0.010000000"] + [f"{t:.9f},1.000000000" for t in range(0, 61)]
    return "\r\n".join(rows) + "\r\n"


def main():
    utf8_console()
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    cfg = load_study(args.config)
    man = pd.read_csv(cfg["run_dir"] / "raw_manifest.csv", encoding="utf-8-sig")
    lockmass = 554.262022 if cfg["polarity"] == "NEG" else 556.276575   # leucine enkephalin [M-H]- / [M+H]+
    bak = PROJECT / f"archive_lmgt_written_by_mzmine_{date.today().isoformat()}"
    n = {"identity": 0, "written": 0, "replaced": 0, "missing": 0, "non_identity": 0}
    for raw in man.raw_path:
        d = Path(raw)
        p = d / "lmgt.inf"
        if not p.exists():
            n["missing"] += 1
            if not args.check:
                p.write_text(identity_text(lockmass), encoding="ascii", newline="")
                n["written"] += 1
        elif is_identity(p):
            n["identity"] += 1
        else:
            n["non_identity"] += 1
            if not args.check:
                bak.mkdir(exist_ok=True)
                shutil.copy2(p, bak / f"{d.name}.lmgt.inf")
                p.write_text(identity_text(lockmass), encoding="ascii", newline="")
                n["replaced"] += 1
    print(("[점검] " if args.check else "[완료] ") + " · ".join(f"{k} {v}" for k, v in n.items()))


if __name__ == "__main__":
    main()
