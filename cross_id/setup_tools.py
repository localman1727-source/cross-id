# -*- coding: utf-8 -*-
"""
새 PC 에 분석 도구 세우기 — MS-DIAL · mzmine
==============================================
'왜' 설치 안내서를 따라 하면 안 되나?
  안내서대로 깔면 **실행이 죽는다.** MS-DIAL 배포본에는 `RDAM_DLL.dll` 과
  `RfxCommonBase.dll` 이 빠져 있어서, 그대로 돌리면 "RDAM_DLL 을 로드할 수 없습니다"로
  끝난다. 이름만 보고는 무엇이 빠졌는지 알기 어렵다. 2026-09-08 에 이 PC 에서 이 문제를
  찾아 해결했고, 그 절차를 그대로 다시 밟는 것이 이 스크립트다.

이 스크립트가 하는 일
  ① 배포본 압축을 풀어 C:\\MSDIAL\\{gui, console} 구조를 만든다 (화면판·콘솔판 둘 다)
  ② 빠진 파일 두 개를 ABF 변환기 폴더에서 찾아 gui·console **양쪽에** 복사한다
  ③ Waters 원본을 변환 없이 읽는 데 필요한 파일이 있는지 본다
  ④ 콘솔을 실제로 한 번 불러 응답하는지 확인한다
  ⑤ mzmine 설치와 로그인 파일을 확인한다 (로그인은 사람이 해야 한다)

하지 않는 일
  인터넷에서 배포본을 내려받지 않는다. 받는 곳이 바뀌거나 판이 올라가면 조용히 다른 것을
  깔게 되기 때문이다. 배포본 압축 파일은 사람이 준비해 경로로 넘긴다.

실행
  py scripts\\setup_tools.py --check                     지금 상태만 본다
  py scripts\\setup_tools.py --msdial-zip <폴더 또는 zip> --abf-dir "<ABF Converter 폴더>"
  py scripts\\setup_tools.py --msdial-dir D:\\MSDIAL ...  다른 위치에 깔 때
"""
import sys
import shutil
import zipfile
import argparse
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import utf8_console, PROJ

utf8_console()

DEFAULT_MSDIAL = Path(r"C:\MSDIAL")
DEFAULT_MZMINE = Path(r"C:\mzmine")
NEED_DLL = ("RDAM_DLL.dll", "RfxCommonBase.dll")
WATERS_DLL = ("MassLynxRaw.dll", "MassLynxSDK.dll")

OK, NO, WARN, DO = "  [있음] ", "  [없음] ", "  [주의] ", "  [작업] "


def find_zips(where: Path):
    """배포본 압축 두 개(화면판 · 콘솔판)를 가려낸다. 파일명에 console 이 있으면 콘솔판이다."""
    if where.is_file():
        cands = [where]
    else:
        cands = sorted(where.glob("MSDIAL*.zip")) or sorted(where.glob("*.zip"))
    gui = next((p for p in cands if "console" not in p.name.lower()), None)
    console = next((p for p in cands if "console" in p.name.lower()), None)
    return gui, console


def unzip_into(zip_path: Path, target: Path) -> bool:
    """압축 안에 폴더가 한 겹 더 있는 경우를 펴서 target 바로 아래에 둔다."""
    if zip_path is None:
        return False
    target.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        tops = {n.split("/")[0] for n in names if "/" in n}
        single_root = len(tops) == 1 and not any("/" not in n and n for n in names[:5])
        z.extractall(target if not single_root else target.parent / "_unzip_tmp")
    if single_root:
        tmp = target.parent / "_unzip_tmp"
        inner = next(tmp.iterdir())
        for item in inner.iterdir():
            dst = target / item.name
            if dst.exists():
                continue
            shutil.move(str(item), str(dst))
        shutil.rmtree(tmp, ignore_errors=True)
    return True


def install_msdial(msdial_dir: Path, zip_src: Path | None, abf_dir: Path | None, check_only: bool):
    print("■ MS-DIAL")
    gui_dir, console_dir = msdial_dir / "gui", msdial_dir / "console"

    if not check_only and zip_src:
        gui_zip, console_zip = find_zips(zip_src)
        if not (gui_zip or console_zip):
            print(f"{NO}배포본 압축을 찾지 못했습니다: {zip_src}")
        for z, d, label in ((gui_zip, gui_dir, "화면판"), (console_zip, console_dir, "콘솔판")):
            exe = d / ("MSDIAL.exe" if label == "화면판" else "MSDIALCUI.exe")
            if exe.exists():
                print(f"{OK}{label} 이미 설치됨  {d}")
            elif z:
                print(f"{DO}{label} 압축 풀기  {z.name} → {d}")
                unzip_into(z, d)
            else:
                print(f"{NO}{label} 압축 파일이 없습니다")

    cui = console_dir / "MSDIALCUI.exe"
    gui = gui_dir / "MSDIAL.exe"
    print(f"{OK if cui.exists() else NO}콘솔  {cui}")
    print(f"{OK if gui.exists() else NO}화면판  {gui}"
          + ("" if gui.exists() else "   ← 파라미터를 내보내려면 필요합니다"))

    # ② 배포본에 없는 파일 두 개
    missing_targets = [d for d in (console_dir, gui_dir) if d.exists()
                       and not all((d / n).exists() for n in NEED_DLL)]
    if not missing_targets:
        if console_dir.exists():
            print(f"{OK}배포본에 없는 파일 두 개 (RDAM_DLL · RfxCommonBase) — 양쪽 다 있음")
    else:
        src_dir = None
        if abf_dir:
            hits = [p for p in Path(abf_dir).rglob(NEED_DLL[0])]
            src_dir = hits[0].parent if hits else None
        if src_dir is None:
            print(f"{NO}RDAM_DLL.dll · RfxCommonBase.dll 을 찾지 못했습니다")
            print("         꾸러미를 만들 때 --with-abf 를 붙이면 이 두 파일이 tools/ 에 들어옵니다.")
            print("         이미 있는 PC 에서는 워크스페이스의 'ABF Converter' 폴더에 들어 있습니다.")
            print("         직접 알려주려면 --abf-dir \"<폴더>\" 를 쓰세요.")
        elif check_only:
            print(f"{WARN}두 파일이 빠져 있습니다 — --check 없이 다시 돌리면 복사합니다")
        else:
            for d in missing_targets:
                for n in NEED_DLL:
                    if not (d / n).exists():
                        shutil.copy2(src_dir / n, d / n)
                        print(f"{DO}복사  {n} → {d.name}\\")

    # ③ Waters 직접 읽기
    if console_dir.exists():
        have = [n for n in WATERS_DLL if (console_dir / n).exists()]
        if len(have) == len(WATERS_DLL):
            print(f"{OK}Waters 원본 직접 읽기 가능 (변환 단계를 건너뛸 수 있습니다)")
        else:
            print(f"{WARN}Waters 직접 읽기용 파일이 부족합니다 — ABF 로 변환해 쓰게 됩니다")

    # ④ 실제로 불러보기
    if cui.exists():
        try:
            # 도구 화면 글자가 한글 윈도우 기본 인코딩이라 utf-8 로 읽으면 터진다 — 바꿔 읽는다
            r = subprocess.run([str(cui), "lcms", "-h"], capture_output=True, text=True,
                               encoding="utf-8", errors="replace", timeout=120)
            print(f"{OK if r.returncode == 0 else NO}콘솔이 응답합니다 (lcms -h)")
        except (OSError, ValueError, subprocess.TimeoutExpired) as e:
            print(f"{NO}콘솔을 부르지 못했습니다: {e}")

    # 메소드 본보기 위치 안내
    ref = PROJ / "msdial" / "reference"
    if ref.exists():
        print(f"{OK}메소드 형식 정본과 지난 기록  {ref}")
        print("         새 메소드를 쓸 때 Parameter_export.txt 의 키 이름을 그대로 씁니다")
    print(f"{WARN}조각무늬 사전(.msp)은 배포본에 없습니다 — 따로 받아 설정의 msdial_msp 에 적습니다")


def install_mzmine(mzmine_dir: Path, zip_src: Path | None, check_only: bool):
    print("\n■ mzmine")
    exe = mzmine_dir / "mzmine_console.exe"
    if not exe.exists() and zip_src and not check_only:
        z = next((p for p in Path(zip_src).glob("*mzmine*.zip")), None)
        if z:
            print(f"{DO}압축 풀기  {z.name} → {mzmine_dir}")
            unzip_into(z, mzmine_dir)
    print(f"{OK if exe.exists() else NO}콘솔  {exe}")

    users = Path.home() / ".mzmine" / "users"
    found = sorted(users.glob("*.mzuser")) if users.exists() else []
    if found:
        print(f"{OK}로그인 파일  {found[0].name}")
    else:
        print(f"{NO}로그인 파일 — mzmine 창을 한 번 띄워 로그인하세요")
        print("         이건 사람이 해야 합니다. 없으면 콘솔이 아무 일도 안 하고 끝납니다")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="설치하지 않고 상태만 본다")
    ap.add_argument("--msdial-zip", default=None, help="MS-DIAL 배포본 zip 이 있는 폴더(또는 zip 파일)")
    ap.add_argument("--msdial-dir", default=str(DEFAULT_MSDIAL), help="MS-DIAL 을 깔 위치")
    ap.add_argument("--mzmine-zip", default=None, help="mzmine 배포본 zip 이 있는 폴더")
    ap.add_argument("--mzmine-dir", default=str(DEFAULT_MZMINE), help="mzmine 을 깔 위치")
    ap.add_argument("--abf-dir", default=None, help="ABF 변환기 폴더 (빠진 파일 두 개를 여기서 가져온다)")
    args = ap.parse_args()

    # ABF 폴더를 안 주면 스스로 찾는다.
    #   순서: 꾸러미의 tools/ → 워크스페이스의 ABF Converter → 없음
    #   '왜' 꾸러미가 먼저인가? 새 PC 에는 워크스페이스가 없다. 꾸러미만 받은 상황이 기본이다.
    abf = Path(args.abf_dir) if args.abf_dir else None
    if abf is None:
        for guess in (PROJ.parent / "tools", PROJ.parent / "ABF Converter",
                      PROJ.parent / "tools" / "ABF Converter"):
            if guess.exists() and next(iter(guess.rglob(NEED_DLL[0])), None):
                abf = guess
                break

    # 배포본 압축도 마찬가지로 꾸러미의 tools/ 를 먼저 본다
    msdial_zip = args.msdial_zip
    if msdial_zip is None:
        t = PROJ.parent / "tools"
        if t.exists() and any(t.glob("MSDIAL*.zip")):
            msdial_zip = str(t)
            print(f"  [안내] 꾸러미 안의 배포본을 씁니다: {t}\n")

    print("=" * 70)
    print(" 분석 도구 구축 — MS-DIAL · mzmine")
    print(f" 자세한 배경은 {PROJ / 'msdial' / 'MSDIAL_구축_규격.md'}")
    print("=" * 70 + "\n")

    install_msdial(Path(args.msdial_dir),
                   Path(msdial_zip) if msdial_zip else None,
                   abf, args.check)
    install_mzmine(Path(args.mzmine_dir),
                   Path(args.mzmine_zip) if args.mzmine_zip else None,
                   args.check)

    print("\n" + "=" * 70)
    print(" 다음: py scripts\\install_assets.py 로 파이썬 쪽까지 한 번에 점검하세요.")
    print("=" * 70)


if __name__ == "__main__":
    main()
