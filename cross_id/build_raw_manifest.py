# -*- coding: utf-8 -*-
"""
① 원시자료 목록 만들기 (어느 연구든 공통)
============================================
'왜' 목록을 파일로 남기나?
  원시자료는 워크스페이스 밖 로컬 디스크에 있고 동기화되지 않는다(규정 ③).
  다른 PC·다음 세션은 "무엇이 몇 개 있었는가"를 이 표로만 알 수 있다.

'왜' 기기 헤더를 읽나?
  파일명은 사람이 붙인 것이라 틀릴 수 있다. 취득 일시·장비명은 기기가 쓴 값이다.
  둘이 어긋나면 그 자체가 경고 신호다.

'왜' 신호 갈래(함수) 수를 세나?
  조각을 한꺼번에 깨는 방식은 보통 세 갈래(전구체 · 조각 · 눈금 맞춤용 기준물질)로 저장된다.
  갈래 수가 예상과 다르면 조각 정보를 못 쓰거나 질량 눈금을 믿을 수 없다는 뜻이다.

산출물: runs/<study_id>/raw_manifest.csv
실행:   py scripts/build_raw_manifest.py --config config/<연구>.json
"""
import sys, csv, re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import load_study, utf8_console

utf8_console()

# 제조사별 원시자료 생김새 — 새 기기를 만나면 여기 한 줄만 는다
VENDOR = {
    "Waters":  {"pattern": "*.raw", "is_dir": True,  "header": "_HEADER.TXT", "func_glob": "_FUNC*.DAT"},
    "Thermo":  {"pattern": "*.raw", "is_dir": False, "header": None,          "func_glob": None},
    "Agilent": {"pattern": "*.d",   "is_dir": True,  "header": None,          "func_glob": None},
    "Shimadzu": {"pattern": "*.qgd", "is_dir": False, "header": None,         "func_glob": None},
}


def read_waters_header(raw_dir: Path) -> dict:
    """'$$ 키: 값' 줄을 사전으로. 기기가 쓴 지역문자가 섞여 있어 errors='replace'."""
    h, f = {}, raw_dir / "_HEADER.TXT"
    if not f.exists():
        return h
    for line in f.read_text(encoding='utf-8', errors='replace').splitlines():
        m = re.match(r"\$\$\s*([^:]+):\s*(.*)", line)
        if m:
            h[m.group(1).strip()] = m.group(2).strip()
    return h


def size_mb(p: Path) -> float:
    if p.is_dir():
        return sum(f.stat().st_size for f in p.rglob('*') if f.is_file()) / 1024 / 1024
    return p.stat().st_size / 1024 / 1024


def main():
    cfg = load_study()
    spec = VENDOR.get(cfg.get("vendor", "Waters"))
    if spec is None:
        raise SystemExit(f"[중단] 모르는 기기: {cfg.get('vendor')} — VENDOR 표에 한 줄 추가하세요")

    root = Path(cfg["raw_dir"])
    if not root.exists():
        raise SystemExit(f"[중단] 원시자료 폴더가 없습니다: {root}\n"
                         f"        다른 PC 라면 그 디스크가 없는 것이 정상입니다 — 산출물만 보고 진행하세요")

    codes = [s["sample_code"] for s in cfg["samples"]]
    rows = []
    for raw in sorted(root.glob(spec["pattern"])):
        name = raw.stem
        h = read_waters_header(raw) if spec["header"] else {}
        n_func = len(list(raw.glob(spec["func_glob"]))) if spec["func_glob"] and raw.is_dir() else ""
        # 파일명 안에서 설정에 적힌 시료 표기를 찾는다. 없으면 반복확인용으로 본다.
        hit = [c for c in codes if c in name]
        rows.append({
            "run_name": name,
            "raw_path": str(raw),
            "sample_code": hit[0] if hit else "QC",
            "is_qc": "N" if hit else "Y",
            "polarity": cfg["polarity"] if f"_{cfg['polarity']}_" in name.upper() else "CHECK",
            "injection": name.split("_")[-1],
            "n_functions": n_func,
            "instrument": h.get("Instrument", ""),
            "acquired_date": h.get("Acquired Date", ""),
            "acquired_time": h.get("Acquired Time", ""),
            "size_MB": round(size_mb(raw), 1),
        })

    if not rows:
        raise SystemExit(f"[중단] {root} 에서 {spec['pattern']} 를 찾지 못했습니다")

    out = cfg["run_dir"] / "raw_manifest.csv"
    with out.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    print(f"[완료] {len(rows)} 런 → {out}")
    n_qc = sum(1 for r in rows if r["is_qc"] == "Y")
    print(f"  시료 {len(rows) - n_qc} 런 · 반복확인용 {n_qc} 런")
    odd = {r["n_functions"] for r in rows if r["n_functions"] != ""}
    if len(odd) > 1:
        print(f"  ⚠ 신호 갈래 수가 런마다 다릅니다 {sorted(odd)} — 같은 방법으로 측정됐는지 확인하세요")
    if any(r["polarity"] == "CHECK" for r in rows):
        print("  ⚠ 파일명에서 극성을 확인하지 못한 런이 있습니다 — 설정의 polarity 와 대조하세요")
    if n_qc == 0:
        print("  ⚠ 반복확인용 시료가 없습니다 — 신호 흔들림을 잴 방법이 없습니다")


if __name__ == "__main__":
    main()
