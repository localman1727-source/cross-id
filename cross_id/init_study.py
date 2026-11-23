# -*- coding: utf-8 -*-
"""
새 연구 설정 만들기 — 원시자료 폴더를 훑어 설정 초안을 쓴다
==============================================================
'왜' 설정을 손으로 쓰지 않게 하나?
  처음 쓰는 사람에게 빈 서식을 주면 어느 칸이 필수인지 모른다. 폴더를 훑어 알 수 있는 것
  (기기 · 극성 · 시료 표기 · 반복 수)은 채워 주고, 사람만 아는 것(학명 · 채집지)은
  「채워주세요」 표시를 남긴다.

🔴 이 스크립트는 추정한다. 추정한 값은 반드시 사람이 확인해야 한다
  파일명 규칙은 실험실마다 다르다. 시료 표기를 잘못 뽑으면 반복 측정이 서로 다른 시료로
  갈리거나 그 반대가 되고, 그 상태로도 분석은 에러 없이 끝까지 돈다.

실행: py scripts/init_study.py "D:\\<원시자료 폴더>" [--name <연구이름>]
"""
import sys
import re
import json
import argparse
from pathlib import Path
from collections import Counter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import utf8_console, CONFIG_DIR

utf8_console()

VENDOR_HINT = [("*.raw", "_FUNC001.DAT", "Waters"), ("*.raw", None, "Thermo"),
               ("*.d", None, "Agilent"), ("*.qgd", None, "Shimadzu")]


def guess_vendor(root: Path):
    # '왜' 첫 파일 하나만 보지 않나? 실패한 런 폴더 하나가 맨 앞에 오면(FUNC 파일이 없는 .raw)
    # 나머지가 멀쩡해도 "형식을 모른다"로 끝났다. 하나라도 조건에 맞으면 그 형식으로 본다.
    for pattern, marker, vendor in VENDOR_HINT:
        hits = list(root.glob(pattern))
        if not hits:
            continue
        if marker:
            if any(h.is_dir() and (h / marker).exists() for h in hits):
                return vendor, pattern
        else:
            if any(h.is_file() for h in hits) or vendor in ("Agilent",):
                return vendor, pattern
    return None, None


def find_raw_folders(root: Path, depth: int = 3):
    """원시자료가 실제로 들어 있는 하위 폴더를 찾는다.

    '왜' 필요한가? 사람들은 보통 상위 폴더(날짜 폴더·프로젝트 폴더)를 알려준다. 그런데 이후 단계는
    전부 「원시자료가 바로 들어 있는 폴더」를 기준으로 파일을 센다. 상위 폴더를 주면 아무것도 못 찾고
    멈췄다(2026-09-30, 다른 PC 에서 실제로 막힘).
    """
    found = []
    for d in [root] + [p for p in root.rglob("*") if p.is_dir()
                        and len(p.relative_to(root).parts) <= depth
                        and p.suffix.lower() not in (".raw", ".d")]:   # .raw/.d 폴더 자체는 런이다
        v, pat = guess_vendor(d)
        if v:
            found.append((d, v, pat))
    return found


def describe_contents(root: Path, depth: int = 3) -> str:
    """못 찾았을 때 무엇이 있었는지 보여준다 — "없다"만 말하면 사람이 다음 행동을 못 정한다."""
    exts = Counter()
    for p in root.rglob("*"):
        if len(p.relative_to(root).parts) > depth:
            continue
        exts[p.suffix.lower() or ("(폴더)" if p.is_dir() else "(확장자 없음)")] += 1
    return ", ".join(f"{e} {n}개" for e, n in exts.most_common(10)) or "(비어 있음)"


def guess_polarity(names):
    up = " ".join(names).upper()
    pos = len(re.findall(r"[_\- ]POS|ES\+|\bP\b", up))
    neg = len(re.findall(r"[_\- ]NEG|ES-", up))
    if pos and not neg:
        return "POS"
    if neg and not pos:
        return "NEG"
    return "CHECK"          # 섞여 있으면 사람이 확인한다 — 반으로 갈라 두 연구로 돌리는 것이 보통이다


def guess_samples(names):
    """파일명에서 시료 표기를 뽑는다.

    규칙: 이름을 밑줄로 쪼개 날짜·극성·주입량·반복번호처럼 보이는 조각을 버리고,
    남은 조각을 이어 붙인 것을 시료 표기로 본다. 여러 런에 같은 표기가 반복되면
    그것이 시료이고, 반복 수는 그 개수다.
    """
    drop = re.compile(r"^(\d{6,8}|POS|NEG|ES[+-]|\d+ul|\d+uL|\d{1,3}|QC|BLANK|STD)$", re.I)
    table = {}
    for n in names:
        parts = [p for p in n.split("_") if p]
        keep = [p for p in parts if not drop.match(p)]
        code = "_".join(keep) if keep else "UNKNOWN"
        if re.search(r"\bQC\b|_QC_|^QC", n, re.I):
            code = "QC"
        table.setdefault(code, []).append(n)
    return table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_dir", help="원시자료가 들어 있는 폴더")
    ap.add_argument("--name", default=None, help="연구 이름 (생략하면 폴더 이름에서 만든다)")
    args = ap.parse_args()

    root = Path(args.raw_dir)
    if not root.exists():
        raise SystemExit(f"[중단] 폴더가 없습니다: {root}")

    vendor, pattern = guess_vendor(root)
    if vendor is None:
        found = find_raw_folders(root)
        if len(found) == 1:
            # 원시자료가 든 하위 폴더가 딱 하나면 그 폴더로 내려간다
            root, vendor, pattern = found[0]
            print(f"[안내] 알려주신 폴더 바로 아래가 아니라 하위 폴더에 원시자료가 있어 그 폴더를 씁니다\n"
                  f"       → {root}\n")
        elif len(found) > 1:
            print(f"[중단] 원시자료가 든 폴더가 {len(found)} 곳입니다. 한 번에 한 폴더(한 연구)만 설정합니다.")
            print("       아래 중 하나를 골라 그 경로로 다시 실행하세요 (보통 양이온·음이온을 따로 잡습니다)\n")
            for d, v, pat in found:
                print(f"   {d}   ({v}, {len(list(d.glob(pat)))} 런)")
            raise SystemExit(1)
        else:
            raise SystemExit(
                f"[중단] {root} 와 그 아래 세 단계까지 뒤졌지만 아는 형식의 원시자료가 없습니다.\n"
                f"        찾는 것: .raw(Waters 폴더/Thermo 파일) · .d(Agilent) · .qgd(Shimadzu)\n"
                f"        실제로 있던 것: {describe_contents(root)}\n"
                f"        .egg · .zip 이 보이면 압축을 먼저 풀고, .mzML 만 있으면 이 절차의 입력이 아닙니다.")

    runs = sorted(p.stem for p in root.glob(pattern))
    polarity = guess_polarity(runs)
    table = guess_samples(runs)
    qc_runs = table.get("QC", [])
    samples = {k: v for k, v in table.items() if k != "QC"}

    study_id = args.name or (re.sub(r"[^0-9A-Za-z가-힣]+", "_", root.name).strip("_") + f"_{polarity}")

    cfg = {
        "study_id": study_id,
        "title": f"{root.name} — 제목을 적어주세요",
        "raw_dir": str(root),
        "vendor": vendor,
        "polarity": polarity,
        "acquisition": "MSe" if vendor == "Waters" else "CHECK — MSe / DDA / DIA 중 무엇인지 적어주세요",
        "samples": [
            {
                "sample_code": code,
                "collect_id": "",
                "name_kr": "채워주세요 — 표준 국명",
                "scientific_name": "채워주세요 — 속 종 (명명자)",
                "genus": "채워주세요 — 속",
                "family": "채워주세요 — 과",
                "note": f"런 {len(v)} 개",
            } for code, v in sorted(samples.items())
        ],
        "family_aliases": ["채워주세요 — 현행 과", "옛 분류의 과 (없으면 지우세요)"],
        "far_genera": [],
        "site_source": {
            "path": "",
            "header_row": 1,
            "columns": {},
            "known_errors": "",
        },
        "processing": {
            "rt_min": 0.3,
            "rt_max": 20,
            "mz_min": 100,
            "mz_max": 1200,
            "min_peak_height": 2000,
            "msdial_msp": "",
            # 극성에 맞는 틀. 전에는 음이온 자료에도 양이온 틀을 박았다 — 에러 없이 틀린 설정으로 돈다
            "mzmine_template": "templates/"
                               + ("example_waters_mse_NEG.mzbatch" if polarity == "NEG" else "example_waters_mse_POS.mzbatch"),
        },
        "design_caveats": [],
    }

    out = CONFIG_DIR / f"{study_id}.json"
    if out.exists():
        raise SystemExit(f"[중단] 같은 이름의 설정이 이미 있습니다: {out}\n"
                         f"        덮어쓰지 않습니다 — --name 으로 다른 이름을 주세요")
    out.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding='utf-8')

    print(f"[생성] {out}\n")
    print(f"  기기       {vendor}")
    print(f"  극성       {polarity}" + ("   ← 확인 필요" if polarity == "CHECK" else ""))
    print(f"  런         모두 {len(runs)} 개 (시료 {len(runs) - len(qc_runs)} · 반복확인용 {len(qc_runs)})")
    print(f"  시료 표기  {', '.join(sorted(samples)) or '(없음)'}")
    rep = Counter(len(v) for v in samples.values())
    print(f"  반복 수    {dict(rep)}")
    if not qc_runs:
        print("\n  ⚠ 반복확인용(QC) 시료를 찾지 못했습니다 — 신호 흔들림을 잴 수 없습니다")
    if len(set(len(v) for v in samples.values())) > 1:
        print("  ⚠ 시료마다 반복 수가 다릅니다 — 파일명 규칙이 섞였는지 확인하세요")

    print("\n다음에 할 일")
    print(f"  1) {out.name} 을 열어 「채워주세요」 칸을 채웁니다 (학명 · 속 · 과)")
    print("  2) 채집지·효능값이 다른 표에 있으면 site_source 를 채웁니다")
    print("  3) 아는 설계상 약점을 design_caveats 에 적습니다 (비우면 보고서에서 사라집니다)")
    print(f"  4) py scripts\\run_stage1.py --config config\\{out.name}")


if __name__ == "__main__":
    main()
