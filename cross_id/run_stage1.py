# -*- coding: utf-8 -*-
"""
1단계 자동 실행 — 원시자료에서 동정 표까지 한 번에
====================================================
무엇을 하나 (그림의 ① ~ ⑥)
  ① 원시자료 목록 · 처리용 시료표          ② 채집지 · 효능값 (있을 때)
  ③ 후보 성분 목록                          ④ mzmine 배치 만들고 돌리기
  ⑤ MS-DIAL 입력 · 메소드 만들고 돌리기     ⑥ 두 결과 맞대보기

'왜' 한 번에 도는 실행기가 필요한가?
  단계가 여덟이고 각 단계가 앞 단계의 산출물을 읽는다. 사람이 순서를 외워 돌리면 한 단계를
  건너뛰고도 다음 단계가 옛 파일을 읽어 에러 없이 돈다 — 결과는 조용히 틀린다.

🔴 오래 걸리는 두 단계는 물어보고 돈다
  두 도구는 런 하나에 수십 분이 걸린다. 처음에는 `--pilot` 으로 시료마다 한 번씩만 돌려
  설정 오류를 털고, 그다음에 전체를 돌린다.

실행
  py scripts\\run_stage1.py --config config\\<연구>.json            준비만 (도구는 안 돌림)
  py scripts\\run_stage1.py --config config\\<연구>.json --run      도구까지 전부
  py scripts\\run_stage1.py --config config\\<연구>.json --run --pilot   시료당 1런만
  py scripts\\run_stage1.py --config config\\<연구>.json --from msdial   중간부터 이어서
"""
import sys
import time
import shutil
import argparse
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study, utf8_console

utf8_console()

MZMINE_EXE = Path(r"C:\mzmine\mzmine_console.exe")
MZUSER = Path.home() / ".mzmine" / "users" / "user.mzuser"
MSDIAL_EXE = Path(r"C:\MSDIAL\console\MSDIALCUI.exe")


def msdial_exe(cfg):
    """설정에 processing.msdial_exe 가 있으면 그 실행 파일을 쓴다.
    '왜'? 정식판 5.5.260820 은 Waters MSe(충돌 세기 여러 개) 자료에서 정렬 뒤 멈추고 런별 결과의
    충돌 세기 선택도 불안정하다(2026-10-08). 고친 자체 빌드를 연구별로 고를 수 있게 둔다."""
    v = cfg["processing"].get("msdial_exe")
    return Path(v) if v else MSDIAL_EXE

STEPS = ["manifest", "metadata", "sites", "candidates", "mzmine", "msdial", "merge", "diagram"]


def say(msg=""):
    print(msg, flush=True)


def run(cmd, log: Path, label: str) -> bool:
    """한 단계를 돌리고 화면에는 요약만, 전문은 기록 파일에 남긴다."""
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with log.open("w", encoding='utf-8', errors='replace') as fh:
        fh.write(" ".join(str(c) for c in cmd) + "\n\n")
        fh.flush()
        p = subprocess.run([str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    ok = p.returncode == 0
    say(f"   {'완료' if ok else '실패'} — {label}  ({dt/60:.1f} 분, 기록 {log.name})")
    if not ok:
        tail = log.read_text(encoding='utf-8', errors='replace').splitlines()[-12:]
        say("   ┌ 마지막 화면")
        for t in tail:
            say(f"   │ {t}")
        say("   └")
    return ok


def preflight(cfg, want_tools: bool) -> bool:
    """돌리기 전에 막는다 — 여기서 막아야 몇 시간짜리 실행이 헛돌지 않는다."""
    say("■ 사전 점검")
    ok = True

    raw = Path(cfg["raw_dir"])
    if raw.exists():
        say(f"   있음   원시자료 폴더 {raw}")
    else:
        say(f"   없음   원시자료 폴더 {raw}  ← 이 PC 에 그 디스크가 있는지 확인하세요")
        ok = False

    todo = [s["sample_code"] for s in cfg["samples"]
            if "채워주세요" in str(s.get("scientific_name", "")) or not s.get("genus")]
    if todo:
        say(f"   미완   설정에 아직 안 채운 칸이 있습니다: {', '.join(todo)}")
        say("          학명·속·과를 채우지 않으면 후보 성분 목록이 비어 나옵니다")
        ok = False
    else:
        say("   있음   시료 설명 (학명 · 속 · 과)")

    tpl = cfg["ws_dir"] / cfg["processing"]["mzmine_template"]
    say(f"   {'있음' if tpl.exists() else '없음'}   mzmine 배치 틀 {tpl.name}")
    ok &= tpl.exists()

    if want_tools:
        for exe, name in ((MZMINE_EXE, "mzmine"), (msdial_exe(cfg), "MS-DIAL")):
            say(f"   {'있음' if exe.exists() else '없음'}   {name} 실행 파일 {exe}")
            ok &= exe.exists()
        say(f"   {'있음' if MZUSER.exists() else '없음'}   mzmine 로그인 파일 {MZUSER.name}"
            + ("" if MZUSER.exists() else "  ← mzmine 을 한 번 띄워 로그인하세요"))
        ok &= MZUSER.exists()

        msp = cfg["processing"].get("msdial_msp", "")
        if msp and Path(msp).exists():
            say(f"   있음   조각무늬 사전 {Path(msp).name}")
        else:
            say("   없음   조각무늬 사전 (MS/MS spectral library, .msp)")
            say("          MS-DIAL 은 돌지만 이름이 하나도 안 붙습니다 — 설정의 msdial_msp 를 채우세요")

        free = shutil.disk_usage("C:\\").free / 1024**3
        say(f"   여유   C 드라이브 {free:.0f} GB" + ("" if free > 50 else "  ← 도구 산출물이 큽니다"))

    say()
    return ok


def write_msdial_method(cfg) -> Path:
    """설정 값으로 콘솔 메소드를 찍어낸다.

    '왜' 틀을 그대로 두고 값만 바꾸나? 콘솔은 모르는 키를 만나면 조용히 무시하고 기본값으로
    돌아간다. 검증된 틀의 줄 구조를 유지하고 값만 갈아끼우는 것이 가장 안전하다.
    """
    tpl = cfg["proj_dir"] / "msdial" / f"method_template_{cfg['polarity']}.txt"
    if not tpl.exists():
        tpl = cfg["proj_dir"] / "msdial" / "method_template_POS.txt"
    pr = cfg["processing"]
    repl = {
        "Ion mode": "Positive" if cfg["polarity"] == "POS" else "Negative",
        "Retention time begin": pr["rt_min"],
        "Retention time end": pr["rt_max"],
        "MS1 mass range begin": pr["mz_min"],
        "MS1 mass range end": pr["mz_max"],
        "Minimum peak height": pr["min_peak_height"],
        "Msp file path": pr.get("msdial_msp", ""),
        "Searched adduct ions": "[M+H]+,[M+Na]+,[M+NH4]+" if cfg["polarity"] == "POS" else "[M-H]-,[M+HCOO]-",
    }
    # 일꾼 수는 설정에 있을 때만 바꾼다 — 사전이 크면 일꾼마다 메모리를 크게 잡는다(2026-10-07 메모리 부족 중단)
    if pr.get("msdial_threads"):
        repl["Number of threads"] = pr["msdial_threads"]
    out_lines = []
    for line in tpl.read_text(encoding='utf-8').splitlines():
        key = line.split(":", 1)[0].strip()
        out_lines.append(f"{key}: {repl[key]}" if key in repl else line)
    out = cfg["run_dir"] / "msdial_method.txt"
    out.write_text("\n".join(out_lines) + "\n", encoding='utf-8')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--run", action="store_true", help="오래 걸리는 두 도구까지 실제로 돌린다")
    ap.add_argument("--pilot", action="store_true", help="시료마다 1런만 (설정 오류 털기용)")
    ap.add_argument("--from", dest="start", default=None, choices=STEPS, help="이 단계부터 이어서")
    ap.add_argument("--only", default=None, choices=STEPS, help="이 단계만")
    args = ap.parse_args()

    cfg = load_study(args.config)
    run_dir, logs = cfg["run_dir"], cfg["run_dir"] / "logs"
    py = sys.executable
    conf = ["--config", str(cfg["config_path"])]

    say("=" * 74)
    say(f" 천연물 비교 연구 1단계 — {cfg.get('title', cfg['study_id'])}")
    say(f" 설정 {cfg['config_path'].name} · 산출물 runs\\{cfg['study_id']}\\")
    say("=" * 74)
    say()

    if not preflight(cfg, want_tools=args.run) and not args.only:
        say("■ 사전 점검에서 막혔습니다. 위 「없음 · 미완」 줄을 해결한 뒤 다시 실행하세요.")
        return 1

    todo = [args.only] if args.only else STEPS[STEPS.index(args.start):] if args.start else STEPS
    done, failed = [], []

    def step(name, label, cmd, skip_reason=None):
        if name not in todo:
            return
        say(f"■ {label}")
        if skip_reason:
            say(f"   건너뜀 — {skip_reason}")
            return
        (done if run(cmd, logs / f"{name}.log", label) else failed).append(name)
        say()

    step("manifest", "① 원시자료 목록", [py, HERE / "build_raw_manifest.py", *conf])
    step("metadata", "① 처리용 시료표", [py, HERE / "build_sample_metadata.py", *conf])
    step("sites", "② 채집지 · 효능값",
         [py, HERE / "build_sample_sites.py", *conf],
         skip_reason=None if any(s.get("path") for s in (cfg.get("site_source") if isinstance(cfg.get("site_source"), list) else [cfg.get("site_source") or {}])) else "설정에 시료 목록이 없습니다")
    step("candidates", "③ 후보 성분 목록", [py, HERE / "build_candidate_db.py", *conf])

    # ── ④ mzmine ───────────────────────────────────────────────────────────
    if "mzmine" in todo:
        say("■ ④ mzmine — 배치 만들기")
        cmd = [py, HERE / "make_mzmine_batches.py", *conf]
        if args.pilot:
            cmd += ["--limit-rep", "1"]
        ok = run(cmd, logs / "mzmine_batch.log", "배치 만들기")
        say()
        batches = sorted((run_dir / "batches").glob("*.mzbatch"))
        if ok and args.run and batches:
            say("■ ④ mzmine — 분석 돌리기 (오래 걸립니다)")
            ok = run([MZMINE_EXE, "-b", batches[-1], "-user", MZUSER, "-m", "all"],
                     logs / "mzmine_run.log", "mzmine 분석")
            (done if ok else failed).append("mzmine")
            say()
        elif ok and not args.run:
            say(f"   준비됨 — 돌리려면 --run 을 붙이세요 (배치 {batches[-1].name if batches else '없음'})")
            say()

    # ── ⑤ MS-DIAL ──────────────────────────────────────────────────────────
    if "msdial" in todo:
        say("■ ⑤ MS-DIAL — 입력 목록과 메소드 만들기")
        ok = run([py, HERE / "make_msdial_input.py", *conf], logs / "msdial_input.log", "입력 목록")
        method = write_msdial_method(cfg)
        say(f"   완료 — 메소드 {method.name}")
        say()
        if ok and args.run:
            out = run_dir / "msdial_out"
            out.mkdir(parents=True, exist_ok=True)
            n = sum(1 for _ in (run_dir / "msdial_input.csv").open(encoding='utf-8')) - 1
            say(f"■ ⑤ MS-DIAL — 분석 돌리기 ({n} 런, 런당 수십 분까지 걸립니다)")
            ok = run([msdial_exe(cfg), "lcms", "-i", run_dir / "msdial_input.csv",
                      "-o", out, "-m", method, "-p"], logs / "msdial_run.log", "MS-DIAL 분석")
            (done if ok else failed).append("msdial")
            say()
        elif ok:
            say("   준비됨 — 돌리려면 --run 을 붙이세요")
            say()

    have_tool_output = (any((run_dir / "mzmine_out").glob("*_full_feature_table.csv"))
                        if (run_dir / "mzmine_out").exists() else False) or                        (any((run_dir / "msdial_out").glob("*.mzTab"))
                        if (run_dir / "msdial_out").exists() else False)
    step("merge", "⑥ 두 결과 맞대보기", [py, HERE / "merge_annotations.py", *conf],
         skip_reason=None if have_tool_output else "두 도구를 아직 안 돌렸습니다 (--run)")
    step("diagram", "모식도", [py, HERE / "make_workflow_diagram.py", *conf])

    # ── 마무리 ─────────────────────────────────────────────────────────────
    say("=" * 74)
    say(f" 끝났습니다. 산출물: runs\\{cfg['study_id']}\\")
    if failed:
        say(f" 실패한 단계: {', '.join(failed)} — logs\\ 의 같은 이름 기록을 보세요")
    if not args.run:
        say(" 도구는 아직 안 돌렸습니다. 준비된 입력으로 돌리려면 --run 을 붙이세요.")
        say(" 처음이라면 --run --pilot 로 시료당 한 런만 먼저 돌려 설정 오류를 터세요.")
    msp = cfg["processing"].get("msdial_msp", "")
    if not (msp and Path(msp).exists()):
        say(" 조각무늬 사전이 없어 MS-DIAL 은 이름을 붙이지 못합니다 (설정의 msdial_msp).")
    say(" 결과를 인용하기 전에 설정의 design_caveats 를 다시 읽으세요 — 이 자료로 말할 수 있는")
    say(" 범위가 거기 적혀 있습니다.")
    say("=" * 74)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
