# -*- coding: utf-8 -*-
"""
연구 설정 읽기 — 이 프로젝트의 모든 스크립트가 공유하는 한 곳
================================================================
'왜' 설정 파일로 빼는가?
  종·산지·기기가 코드 안에 박혀 있으면 다음 연구에서 코드를 복사해 고치게 된다.
  그러면 사본이 늘고, 한쪽만 고쳐진 채 남아 조용히 다른 결과를 낸다(워크스페이스 규정 3-0).
  값은 config/ 의 파일에, 절차는 scripts/ 의 코드에 둔다. 코드는 한 벌만 존재한다.

'왜' 경로를 스스로 계산하나?
  계정명이나 저장 위치가 박히면 다른 PC 에서 깨진다(규정 3-1). 이 파일 위치에서
  프로젝트 뿌리와 워크스페이스 뿌리를 거슬러 찾는다.

쓰는 법:
    from study_config import load_study
    cfg = load_study()                       # --config 인자 또는 단일 설정 자동 선택
    cfg["run_dir"] / "raw_manifest.csv"
"""
import sys
import io
import json
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
WS = PROJ                             # repository root — relative paths in the config (e.g. templates/) resolve from here
CONFIG_DIR = PROJ / "config"


def utf8_console():
    """윈도우 콘솔 한글 깨짐 방지 (규정 ②). 스크립트 첫 줄에서 부른다."""
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')


def add_config_arg(ap: argparse.ArgumentParser):
    ap.add_argument("--config", default=None,
                    help="config/ 안의 연구 설정 파일. 생략하면 설정이 하나뿐일 때 그것을 쓴다")
    return ap


def load_study(config_path=None) -> dict:
    """설정을 읽고 경로를 붙여 돌려준다.

    설정이 여러 개인데 --config 를 안 주면 멈춘다. '아무거나 골라 도는 것'이
    가장 나쁜 실패다 — 엉뚱한 연구의 산출물을 덮어쓰기 때문이다.
    """
    if config_path is None:
        # 명령줄에서 직접 읽는다 (각 스크립트가 argparse 를 따로 쓰더라도 동작하게)
        for i, a in enumerate(sys.argv):
            if a == "--config" and i + 1 < len(sys.argv):
                config_path = sys.argv[i + 1]
                break
    if config_path is None:
        cands = [p for p in sorted(CONFIG_DIR.glob("*.json"))
                 if p.name != "study_template.json"]
        if len(cands) != 1:
            names = ", ".join(p.name for p in cands) or "(없음)"
            raise SystemExit(
                f"[중단] 연구 설정을 고를 수 없습니다. --config 로 지정하세요.\n"
                f"        config/ 안의 설정: {names}")
        config_path = cands[0]

    p = Path(config_path)
    if not p.is_absolute():
        p = (CONFIG_DIR / p.name) if not p.exists() else p
    if not p.exists():
        raise SystemExit(f"[중단] 설정 파일이 없습니다: {p}")

    cfg = json.loads(p.read_text(encoding='utf-8'))
    cfg = {k: v for k, v in cfg.items() if not k.startswith("_")}   # 설명용 키는 버린다

    for key in ("study_id", "polarity", "samples"):
        if not cfg.get(key):
            raise SystemExit(f"[중단] 설정에 '{key}' 가 비어 있습니다: {p.name}")

    cfg["config_path"] = p
    cfg["proj_dir"] = PROJ
    cfg["ws_dir"] = WS
    cfg["run_dir"] = PROJ / "runs" / cfg["study_id"]
    cfg["db_dir"] = cfg["run_dir"] / "db"
    cfg["run_dir"].mkdir(parents=True, exist_ok=True)
    cfg["db_dir"].mkdir(parents=True, exist_ok=True)
    return cfg


def sample_codes(cfg) -> list:
    return [s["sample_code"] for s in cfg["samples"]]


def taxa_rows(cfg) -> list:
    """(sample_code, collect_id, 국명, 학명, 속, 과, 비고) 형태로 돌려준다."""
    return [(s["sample_code"], s.get("collect_id", ""), s.get("name_kr", ""),
             s.get("scientific_name", ""), s.get("genus", ""), s.get("family", ""),
             s.get("note", "")) for s in cfg["samples"]]
