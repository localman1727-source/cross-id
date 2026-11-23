# -*- coding: utf-8 -*-
"""
SIRIUS 로그인 상태 점검 (2026-11-11, 사무실 큐용)
==============================================
'왜' 따로 있나? SIRIUS 는 로그인이 풀려 있어도 계산을 exit 0 으로 끝낸다(2026-11-06 실사고).
  긴 계산을 걸기 전에 사무실 PC 의 로그인이 살아 있는지만 따로 보고, 결과를 큐 완료 메일로 받는다.
판정: `login --show` 에 'Logged in as' 가 있고 Login ERROR·Not Logged in·invalid_grant 가 없으면 정상.
  `login --limits` 의 라이선스 문구도 함께 찍는다. 아무것도 바꾸지 않는다(읽기만).
종료코드: 0 = 로그인 살아 있음 · 1 = 풀림/오류 · 2 = SIRIUS 없음
쓰는 법: py 천연물_비교연구_project/scripts/check_sirius_login.py
"""
import sys, io, re, subprocess
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

SIRIUS = Path("C:/SIRIUS/sirius/sirius.exe")   # 스킬 natural-product-annotation 0절과 같은 위치


def ask(*args):
    """sirius.exe 하위 명령 하나를 돌려 출력 전체를 돌려준다."""
    r = subprocess.run([str(SIRIUS), *args], capture_output=True, text=True, errors="replace", timeout=300)
    return f"{r.stdout or ''}{r.stderr or ''}"


def main():
    if not SIRIUS.exists():
        print(f"⛔ SIRIUS 없음: {SIRIUS}")
        return 2
    show = ask("login", "--show")
    limits = ask("login", "--limits")
    print("── login --show ──\n" + show.strip()[-2000:])
    print("── login --limits ──\n" + limits.strip()[-1000:])
    bad = re.search(r"Login ERROR|Not Logged in|invalid_grant", show + limits)
    ok = "Logged in as" in show and not bad
    if ok:
        print("✅ SIRIUS 로그인 살아 있음")
        return 0
    print(f"⛔ SIRIUS 로그인 풀림 ({bad.group(0) if bad else 'Logged in as 없음'}) — "
          "C:\\SIRIUS\\sirius\\sirius-gui.exe 에서 로그인 필요")
    return 1


if __name__ == "__main__":
    sys.exit(main())
