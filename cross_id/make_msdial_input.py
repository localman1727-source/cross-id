# -*- coding: utf-8 -*-
"""
④ MS-DIAL 콘솔 입력 목록 만들기
=================================
'왜' 손으로 안 쓰나?
  경로 오타 하나에 그 런만 조용히 빠진 채 끝까지 돈다. 원시자료 목록에서 뽑으면
  실제로 있는 것만 들어간다.

'왜' 반복확인용 시료를 따로 표시하나?
  그 표시가 있어야 도구가 신호 흔들림 보정을 걸 수 있다. 표시를 빼면 그 시료가
  평범한 시료로 섞여 비교 결과를 흐린다.

'왜' 주입 순서를 적나?
  순서를 적어야 나중에 "시간이 지나며 신호가 내려갔는가"를 물을 수 있다.

산출물: runs/<study_id>/msdial_input.csv
실행:   py scripts/make_msdial_input.py --config config/<연구>.json [--abf-dir <변환본 폴더>]
"""
import sys, csv, argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from study_config import load_study, utf8_console

utf8_console()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--abf-dir", default=None,
                    help="변환본(.abf) 폴더. 주면 원본 대신 변환본을 가리킨다")
    args = ap.parse_args()

    cfg = load_study(args.config)
    man = cfg["run_dir"] / "raw_manifest.csv"
    if not man.exists():
        raise SystemExit(f"[중단] 먼저 build_raw_manifest.py 를 돌리세요: {man}")

    rows = list(csv.DictReader(man.open(encoding="utf-8-sig")))
    rows.sort(key=lambda r: (r.get("acquired_date", ""), r.get("acquired_time", "")))

    # 조각을 한꺼번에 깨는 방식은 AIF 로 알려준다. 방식이 다르면 여기서 갈린다.
    acq = {"MSe": "AIF", "DIA": "AIF", "DDA": "DDA", "MS1only": "SCAN"}.get(cfg.get("acquisition", "MSe"), "AIF")

    out_rows = []
    for i, r in enumerate(rows, start=1):
        name = r["run_name"]
        if args.abf_dir:
            path, ftype = str(Path(args.abf_dir) / (name + ".abf")), "ABF"
        else:
            path, ftype = r["raw_path"], {"Waters": "RAW", "Thermo": "RAW",
                                          "Agilent": "ABF", "Shimadzu": "ABF"}.get(cfg.get("vendor"), "RAW")
        out_rows.append([path, name, ftype, r["sample_code"], acq, 1, i, "NA"])

    out = cfg["run_dir"] / "msdial_input.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["file_path", "file_name", "file_type", "class_id",
                    "acquisition_type", "batch_order", "analytical_order", "factor"])
        w.writerows(out_rows)
    print(f"[생성] {out}  ({len(out_rows)} 런 · 방식 {acq})")
    for o in out_rows:
        print(f"  {o[6]:>2}  {o[3]:<10} {o[1]}")


if __name__ == "__main__":
    main()
