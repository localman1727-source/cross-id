# -*- coding: utf-8 -*-
"""
①-b 처리용 시료표 만들기 (원시자료 목록 → 도구가 읽는 표)
============================================================
'왜' 목록을 한 번 더 가공하나?
  원시자료 목록은 "디스크에 무엇이 있나"를 적은 것이고, 이 표는 "그것을 어떻게 묶어
  분석하나"를 적은 것이다. 둘을 한 파일로 합치면, 묶는 방식을 바꿀 때 원본 기록까지
  덮어쓰게 된다.

'왜' 반복 번호를 파일명 끝에서 뽑나?
  같은 시료를 여러 번 주입한 것을 도구가 알아야 흩어진 정도를 잴 수 있다. 이 표에서
  묶어주지 않으면 도구는 서로 다른 시료 아홉 개로 본다.

산출물: runs/<study_id>/sample_metadata.csv
실행:   py scripts/build_sample_metadata.py --config config/<연구>.json
"""
import sys, csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import load_study, utf8_console

utf8_console()


def main():
    cfg = load_study()
    man = cfg["run_dir"] / "raw_manifest.csv"
    if not man.exists():
        raise SystemExit(f"[중단] 먼저 build_raw_manifest.py 를 돌리세요: {man}")

    rows_in = list(csv.DictReader(man.open(encoding="utf-8-sig")))
    name_kr = {s["sample_code"]: s.get("name_kr", "") for s in cfg["samples"]}

    out_rows = []
    for r in rows_in:
        code = r["sample_code"]
        out_rows.append({
            "raw_name": r["run_name"],
            "raw_path": r["raw_path"],
            "sample_code": code,
            "sample_name": code,
            "species_kr": name_kr.get(code, "반복확인용"),
            "sample_note": "",
            "polarity": r["polarity"],
            "replicate": r["injection"].lstrip("0") or "1",
            "acq_date": (r.get("acquired_date") or "").replace("-", ""),
            "analytical_order": "",
        })

    # 실제 주입 순서 = 취득 시각 순. 순서 효과를 보려면 이 값이 맞아야 한다.
    order = sorted(range(len(rows_in)),
                   key=lambda i: (rows_in[i].get("acquired_date", ""), rows_in[i].get("acquired_time", "")))
    for rank, i in enumerate(order, start=1):
        out_rows[i]["analytical_order"] = rank

    out = cfg["run_dir"] / "sample_metadata.csv"
    with out.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)
    print(f"[생성] {out}  ({len(out_rows)} 줄)")

    # 시료가 순서대로 묶여 들어갔는지 — 묶여 있으면 시료 차이와 순서 효과가 겹친다
    seq = [out_rows[i]["sample_code"] for i in order]
    blocks = sum(1 for a, b in zip(seq, seq[1:]) if a != b) + 1
    if blocks <= len(set(seq)):
        print("  ⚠ 시료가 종류별로 묶여 연달아 주입됐습니다 — 시료 차이와 주입 순서가 겹칩니다")
    qc_pos = [k for k, c in enumerate(seq, start=1) if c == "QC"]
    if qc_pos and max(qc_pos) <= len(seq) / 3:
        print("  ⚠ 반복확인용 시료가 앞쪽에만 있습니다 — 시간에 따른 흔들림을 가를 수 없습니다")


if __name__ == "__main__":
    main()
