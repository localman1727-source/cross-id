# -*- coding: utf-8 -*-
"""
학생용 메뉴 — 명령어를 외우지 않고 번호로 고른다
==================================================
'왜' 메뉴를 따로 두나?
  명령줄 인자를 외워야 하면 쓰는 사람이 준다. 그리고 인자를 잘못 주면 조용히 다른 연구의
  산출물을 덮어쓴다. 메뉴는 있는 설정을 보여주고 고르게 하므로 그 사고가 생기지 않는다.

'왜' 메뉴 글씨를 .bat 가 아니라 파이썬이 찍나?
  윈도우 명령창은 배치 파일의 한글을 코드페이지에 따라 깨뜨린다. 파이썬이 찍으면 깨지지 않는다.

실행: 프로젝트 폴더의 「1단계 실행하기.bat」 을 두 번 누르거나  py scripts/menu.py
"""
import sys
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import utf8_console, CONFIG_DIR, PROJ

utf8_console()
PY = sys.executable


def ask(prompt, default=None):
    try:
        v = input(prompt).strip()
    except EOFError:
        return default
    return v or default


def pick_config():
    cands = [p for p in sorted(CONFIG_DIR.glob("*.json")) if p.name != "study_template.json"]
    if not cands:
        print("\n  설정이 하나도 없습니다. 먼저 1번(새 연구 설정 만들기)을 하세요.")
        return None
    print("\n  어느 연구입니까?")
    for i, p in enumerate(cands, 1):
        print(f"   {i}) {p.stem}")
    v = ask("  번호: ")
    if not (v and v.isdigit() and 1 <= int(v) <= len(cands)):
        print("  취소했습니다.")
        return None
    return cands[int(v) - 1]


def call(script, *args):
    cmd = [PY, str(HERE / script), *[str(a) for a in args]]
    print()
    subprocess.run(cmd)
    print()


def main():
    while True:
        print("=" * 66)
        print(" 천연물 비교 연구 — 1단계 실행기")
        print("=" * 66)
        print(" 1) 새 연구 설정 만들기   (원시자료 폴더를 훑어 초안을 씁니다)")
        print(" 2) 준비만 해보기         (목록·후보목록·배치까지. 도구는 안 돌립니다)")
        print(" 3) 시험 삼아 돌리기      (시료마다 1런만. 설정 오류를 먼저 텁니다)")
        print(" 4) 전부 돌리기           (오래 걸립니다. 런당 수십 분)")
        print(" 5) 두 결과 합치기만")
        print(" 6) 작업 모식도 그리기")
        print(" 7) 다른 PC 로 옮길 꾸러미 만들기")
        print(" 8) 이 PC 상태 점검   (무엇이 깔려 있고 무엇이 없는지)")
        print(" 9) 새 PC 에 분석 도구 세우기   (MS-DIAL · mzmine)")
        print(" 0) 끝내기")
        choice = ask("\n 번호를 고르세요: ")

        if choice == "1":
            raw = ask("\n  원시자료 폴더 경로를 넣으세요 (창에 폴더를 끌어다 놓아도 됩니다)\n  > ")
            if raw:
                call("init_study.py", raw.strip('"'))
                print("  설정 파일을 열어 「채워주세요」 칸을 채운 뒤 2번으로 오세요.")
        elif choice in ("2", "3", "4"):
            cfg = pick_config()
            if cfg:
                extra = {"2": [], "3": ["--run", "--pilot"], "4": ["--run"]}[choice]
                if choice == "4":
                    if ask("\n  런 하나에 수십 분씩 걸립니다. 계속할까요? (y/n) ", "n").lower() != "y":
                        continue
                call("run_stage1.py", "--config", cfg, *extra)
        elif choice == "5":
            cfg = pick_config()
            if cfg:
                call("merge_annotations.py", "--config", cfg)
        elif choice == "6":
            cfg = pick_config()
            if cfg:
                call("make_workflow_diagram.py", "--config", cfg)
        elif choice == "7":
            print("\n  원시자료는 너무 커서 꾸러미에 넣지 않습니다 — 외장 디스크로 따로 옮기세요.")
            runs = ask("  지금까지의 결과 표도 넣을까요? (y/n) ", "n").lower() == "y"
            cache = ask("  공개 천연물 자료(약 106 MB)도 넣을까요? 받는 PC 에 인터넷이 없으면 넣습니다 (y/n) ", "n").lower() == "y"
            extra = (["--with-runs"] if runs else []) + (["--with-cache"] if cache else [])
            call("package_assets.py", *extra)
        elif choice == "8":
            call("install_assets.py")
        elif choice == "9":
            print("\n  MS-DIAL 배포본 압축(화면판 · 콘솔판)이 있는 폴더를 넣으세요.")
            print("  그냥 엔터를 치면 지금 상태만 점검합니다.")
            zips = ask("  > ")
            if zips:
                call("setup_tools.py", "--msdial-zip", zips.strip('"'))
            else:
                call("setup_tools.py", "--check")
        elif choice in ("0", "", None):
            print("\n 끝냅니다.")
            return 0
        else:
            print("\n  그런 번호는 없습니다.")

        print(f" 산출물은 {PROJ / 'runs'} 아래에 있습니다.")
        ask(" 계속하려면 엔터를 누르세요. ")
        print()


if __name__ == "__main__":
    raise SystemExit(main())
