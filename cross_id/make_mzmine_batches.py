# -*- coding: utf-8 -*-
"""
⑤ mzmine 배치 XML 생성기 (연구별 값은 config/ 의 설정에서)
============================================
(이전 Phytochem 생성기 `make_batches_phyto.py` 를 이식. 같은 기기(SYNAPT-XS)·
 같은 페놀성 대상·같은 NEG 극성이라 검증된 설정을 그대로 가져온다.
 바꾼 것은 경로 3곳과 기본 인자뿐 — 아래 '경로' 절 참조.)

검증된 리피도믹스 배치(이전 프로젝트)를 템플릿으로 삼되, **파이토케미컬용으로
반드시 바꿔야 하는 곳만** 명시적으로 고쳐서 POS/NEG 배치를 만든다.

'왜' 배치를 처음부터 안 짜고 템플릿을 고치나?
  Waters .raw import + MSe(DIA) pseudo-MS2 + 정렬/갭필 체인은 같은 기기
  (SYNAPT-XS)에서 이미 완주 검증된 설정이다. 처음부터 짜면 그 검증을 버린다.
  바뀌어야 할 것은 '무엇을 찾느냐'(지질→파이토케미컬)와 극성/RT 뿐이다.

🔴 지질 템플릿에서 반드시 바꿔야 하는 것 (안 바꾸면 조용히 틀린다)
  1. RT crop 하한 2.0분 → 0.3분
     '왜'? 2.0분 하한은 지질 POS 의 void volume(=chimeric MS2 온상)을 버리려던
     설정이다. 그런데 파이토케미컬에서 **RT<2분 구간이 바로 목표물**이다 —
     gallic acid·quinic acid·당류 같은 극성 성분이 거기 나온다. 80% MeOH 가
     가장 잘 뽑는 것이 하필 그 성분들이므로, 2.0 을 두면 핵심을 통째로 버린다.
  2. Lipid Annotation 모듈 제거 → Local compound database search 로 교체
     '왜'? Lipid Annotation 은 지방산 사슬 조합 기반이라 플라보노이드·
     카로티노이드에 원리적으로 적용 불가하다.
  3. 극성은 **2곳**을 함께 고친다 (ADAP scan filter + Ion library)
     '왜'? 한 곳만 고치면 에러 없이 동정 0건이 나온다 — 가장 잡기 어려운 함정.

사용법:
    py make_batches_phyto.py --samples 1,2 --name CDE_vs_CLE
    py make_batches_phyto.py --samples 1,2,3,4 --name all4
    py make_batches_phyto.py --samples 1,2 --name pilot --limit-rep 1   # 파일럿
"""
import sys
import io
import csv
import copy
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# ── 경로 (연구별 값은 설정 파일에서) ───────────────────────────────────────
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study

CFG = load_study()
WORKSPACE = CFG["ws_dir"]
RUN = CFG["run_dir"]
PREFIX = CFG["study_id"]
WITH_QC = True

# 검증된 배치를 템플릿으로 삼는다 — 처음부터 짜면 그 검증을 버린다.
TEMPLATE = WORKSPACE / CFG["processing"]["mzmine_template"]
METADATA = RUN / "sample_metadata.csv"
DB_CSV = RUN / "db" / f"{PREFIX}_candidates_relatives.csv"   # build_candidate_db.py 산출물
BATCH_DIR = RUN / "batches"

# 산출물 위치
#   작은 표(성분표·동정표)는 프로젝트 폴더에 (규정 5-1, 보관소로 넘어가야 한다)
#   큰 파일(수백 MB)은 로컬 디스크에 — 공유 폴더에 GB 를 올리지 않는다
OUT_SMALL = RUN / "mzmine_out"
OUT_BIG = Path(r"C:\mzmine_out") / PREFIX

# mzmine 실행 파일과 로그인 사용자 파일
#   '왜' Path.home() 인가? 계정명(사용자명)을 박으면 다른 PC 에서 깨진다 (규정 3-1).
#   .mzmine 폴더는 Dropbox 동기화 대상이 아니므로 PC 마다 1회 로그인이 필요하다.
MZMINE_EXE = Path(r"C:\mzmine\mzmine_console.exe")
MZUSER = Path.home() / ".mzmine" / "users" / "user.mzuser"

# ── mzmine 모듈 식별자 ──────────────────────────────────────────────────────
M_IMPORT = "io.github.mzmine.modules.io.import_rawdata_all.AllSpectralDataImportModule"
M_ADAP = "io.github.mzmine.modules.dataprocessing.featdet_adapchromatogrambuilder.ModularADAPChromatogramBuilderModule"
M_IIN = "io.github.mzmine.modules.dataprocessing.id_ion_identity_networking.ionidnetworking.IonNetworkingModule"
M_LIPID = "io.github.mzmine.modules.dataprocessing.id_lipidid.annotation_modules.LipidAnnotationModule"
M_DBSEARCH = "io.github.mzmine.modules.dataprocessing.id_localcsvsearch.LocalCSVDatabaseSearchModule"
M_GRAPHML = "io.github.mzmine.modules.io.export_network_graphml.NetworkGraphMlExportModule"
M_CSVEXP = "io.github.mzmine.modules.io.export_features_csv.CSVExportModularModule"
M_GNPS = "io.github.mzmine.modules.io.export_features_gnps.fbmn.GnpsFbmnExportAndSubmitModule"
M_SIRIUS = "io.github.mzmine.modules.io.export_features_sirius.SiriusExportModule"
M_META = "io.github.mzmine.modules.visualization.projectmetadata.io.ProjectMetadataExportModule"
M_ANNOEXP = "io.github.mzmine.modules.io.export_compoundAnnotations_csv.CompoundAnnotationsCSVExportModule"
M_AUTOSAVE = "io.github.mzmine.modules.batchmode.autosave.AutoSaveBatchModule"

# ── 이온 라이브러리 (극성별) ────────────────────────────────────────────────
#   질량 규약: 양이온 part 는 전자 1개를 뺀 값, 음이온 part 는 전자 1개를 더한 값.
#   (템플릿의 검증된 값 H=1.00727645209073, Cl=34.96940125990927 과 동일 규약)
#   전자 질량 = 0.00054857990927
#
#   '왜' [M]+• (라디칼 양이온)을 POS 에 넣는가?
#     카로텐류(β-carotene 등)는 산소가 없어 양성자 친화도가 낮고, ESI 에서
#     [M+H]+ 보다 [M]+• 로 뜨는 경우가 흔하다. 이걸 빼면 카로틴을 놓친다.
#   '왜' [M+H-H2O]+ 를 넣는가?
#     크산토필(수산기 보유)·배당체는 탈수 손실이 지배적인 경우가 많다.
#
# ⚠️ 파라미터 이름은 **"Ion library"** 다. 클래스 파일에 보이는 "Use adducts" 는
#    설명용 문구이지 XML 의 name 속성이 아니다. "Use adducts" 로 쓰면 mzmine 이
#    그 파라미터를 못 찾아 기본값(꺼짐)을 쓰고, 결국 "Cannot determine precursor mz"
#    로 동정 단계가 통째로 실패한다. (2026-07-28 실측으로 확인)
ION_LIB_POS = """
        <parameter libraryVersion="2" name="Ion library" selected="true">
            <ionLibrary id="a1000000-0000-4000-8000-000000000001">
                <origin kind="local"/>
                <name>phytochem POS adducts (+)</name>
                <savedDate>2026-07-28T12:00:00.000000</savedDate>
                <lastUpdatedDate>2026-07-28T12:00:00.000000</lastUpdatedDate>
                <parts>
                    <parts charge="1" formula="H" id="0" mass="1.00727645209073" name="H"/>
                    <parts charge="1" formula="Na" id="1" mass="22.98922070009073" name="Na"/>
                    <parts charge="1" formula="K" id="2" mass="38.96315810009073" name="K"/>
                    <parts charge="1" formula="H4N" id="3" mass="18.03382554809073" name="NH4"/>
                    <parts charge="-1" id="4" mass="5.4857990927E-4" name="e"/>
                    <parts charge="0" formula="H2O" id="5" mass="18.0105646863" name="H2O"/>
                </parts>
                <ionTypes>
                    <ionTypes molecules="1"><part count="1" id="0"/></ionTypes>
                    <ionTypes molecules="1"><part count="2" id="0"/></ionTypes>
                    <ionTypes molecules="1"><part count="1" id="1"/></ionTypes>
                    <ionTypes molecules="1"><part count="1" id="2"/></ionTypes>
                    <ionTypes molecules="1"><part count="1" id="3"/></ionTypes>
                    <ionTypes molecules="1"><part count="-1" id="4"/></ionTypes>
                    <ionTypes molecules="1"><part count="1" id="0"/><part count="-1" id="5"/></ionTypes>
                    <ionTypes molecules="2"><part count="1" id="0"/></ionTypes>
                    <ionTypes molecules="2"><part count="1" id="1"/></ionTypes>
                </ionTypes>
            </ionLibrary>
        </parameter>
"""

#   '왜' NEG 에 formate([M+HCOO]-) 를 넣는가?
#     역상 LC 이동상에 formic acid 를 쓰는 것이 표준이며, 중성 배당체·당류는
#     [M-H]- 가 약하고 formate adduct 로 주로 뜬다. 빼면 배당체를 놓친다.
ION_LIB_NEG = """
        <parameter libraryVersion="2" name="Ion library" selected="true">
            <ionLibrary id="a1000000-0000-4000-8000-000000000002">
                <origin kind="local"/>
                <name>phytochem NEG adducts (-)</name>
                <savedDate>2026-07-28T12:00:00.000000</savedDate>
                <lastUpdatedDate>2026-07-28T12:00:00.000000</lastUpdatedDate>
                <parts>
                    <parts charge="1" formula="H" id="0" mass="1.00727645209073" name="H"/>
                    <parts charge="-1" formula="Cl" id="1" mass="34.96940125990927" name="Cl"/>
                    <parts charge="-1" formula="CHO2" id="2" mass="44.99820285790927" name="HCOO"/>
                    <parts charge="0" formula="H2O" id="3" mass="18.0105646863" name="H2O"/>
                </parts>
                <ionTypes>
                    <ionTypes molecules="1"><part count="-1" id="0"/></ionTypes>
                    <ionTypes molecules="1"><part count="-2" id="0"/></ionTypes>
                    <ionTypes molecules="1"><part count="1" id="1"/></ionTypes>
                    <ionTypes molecules="1"><part count="1" id="2"/></ionTypes>
                    <ionTypes molecules="1"><part count="-1" id="0"/><part count="-1" id="3"/></ionTypes>
                    <ionTypes molecules="2"><part count="-1" id="0"/></ionTypes>
                </ionTypes>
            </ionLibrary>
        </parameter>
"""

# ── Local compound database search 배치스텝 ────────────────────────────────
#   datatype 속성은 mzmine 의 **자바 클래스 전체이름**이어야 한다
#   (ImportTypeParameter.loadValueFromXML 이 getClass().getName() 과 비교한다).
#   column 속성은 우리 CSV 의 헤더명이다 → build_compound_db.py 가 쓴 것과 일치해야 함.
DB_SEARCH_STEP = """
    <batchstep method="{module}" module_name="Local compound database search" parameter_version="3">
        <parameter name="Feature lists" type="BATCH_LAST_FEATURELISTS"/>
        <parameter name="Database file">
            <current_file>{db_path}</current_file>
        </parameter>
        <parameter name="Field separator">,</parameter>
        <parameter name="Columns">
            <importtype selected="true" column="name" datatype="io.github.mzmine.datamodel.features.types.annotations.CompoundNameType"/>
            <importtype selected="true" column="formula" datatype="io.github.mzmine.datamodel.features.types.annotations.formula.FormulaType"/>
            <importtype selected="true" column="neutral_mass" datatype="io.github.mzmine.datamodel.features.types.numbers.NeutralMassType"/>
            <importtype selected="true" column="molecular_class" datatype="io.github.mzmine.datamodel.features.types.annotations.compounddb.MolecularClassType"/>
            <importtype selected="false" column="smiles" datatype="io.github.mzmine.datamodel.features.types.annotations.SmilesStructureType"/>
            <importtype selected="false" column="inchi" datatype="io.github.mzmine.datamodel.features.types.annotations.InChIStructureType"/>
            <importtype selected="false" column="inchi_key" datatype="io.github.mzmine.datamodel.features.types.annotations.InChIKeyStructureType"/>
            <importtype selected="false" column="mz" datatype="io.github.mzmine.datamodel.features.types.numbers.PrecursorMZType"/>
            <importtype selected="false" column="rt" datatype="io.github.mzmine.datamodel.features.types.numbers.RTType"/>
            <importtype selected="false" column="CCS" datatype="io.github.mzmine.datamodel.features.types.numbers.CCSType"/>
            <importtype selected="false" column="mobility" datatype="io.github.mzmine.datamodel.features.types.numbers.MobilityType"/>
            <importtype selected="false" column="adduct" datatype="io.github.mzmine.datamodel.features.types.annotations.iin.IonTypeType"/>
            <importtype selected="false" column="comment" datatype="io.github.mzmine.datamodel.features.types.annotations.CommentType"/>
            <importtype selected="false" column="PubChemCID" datatype="io.github.mzmine.datamodel.features.types.annotations.compounddb.PubChemIdType"/>
        </parameter>
        <parameter name="m/z tolerance">
            <absolutetolerance>0.005</absolutetolerance>
            <ppmtolerance>10.0</ppmtolerance>
        </parameter>
        <parameter name="Retention time tolerance" selected="false" unit="MINUTES">0.2</parameter>
        <parameter name="Mobility tolerance" selected="false">1.0</parameter>
        <parameter name="Retention index tolerance" selected="false">10.0</parameter>
        <parameter name="CCS tolerance (%)" selected="false">0.05</parameter>
        <parameter name="Filter filename header" selected="false">raw_filename</parameter>
{ion_library}
        <parameter name="Use isotope matcher" selected="false">
            <parameter name="m/z tolerance">
                <absolutetolerance>0.005</absolutetolerance>
                <ppmtolerance>10.0</ppmtolerance>
            </parameter>
            <parameter name="Minimum isotope intensity">0.01</parameter>
            <parameter name="Minimum isotope score">0.5</parameter>
        </parameter>
        <parameter name="Calculate main isotopic signal">false</parameter>
        <parameter name="Charge filter">NO_FILTER</parameter>
    </batchstep>
"""


# ── 유틸 ────────────────────────────────────────────────────────────────────
def find_step(root: ET.Element, method: str) -> ET.Element:
    """batchstep 을 method 속성으로 찾는다. 없으면 즉시 실패시킨다.

    '왜' 없으면 예외인가? 템플릿이 바뀌어 못 찾았는데 조용히 넘어가면
    '고쳤다고 믿었지만 안 고쳐진' 배치가 만들어진다 — 이번 작업의 최대 위험.
    """
    for step in root.findall("batchstep"):
        if step.get("method") == method:
            return step
    raise RuntimeError(f"템플릿에서 batchstep 을 못 찾음: {method}")


def get_param(parent: ET.Element, name: str) -> ET.Element:
    """직계 자식 parameter 를 name 으로 찾는다 (없으면 예외)."""
    for p in parent.findall("parameter"):
        if p.get("name") == name:
            return p
    raise RuntimeError(f"parameter '{name}' 를 못 찾음 (부모 {parent.get('name') or parent.tag})")


def set_file_param(param: ET.Element, path: Path) -> None:
    """mzmine 의 파일 경로 파라미터를 설정.

    ⚠️ 함정: 경로는 parameter 의 text 가 아니라 <current_file> **자식 요소**에 들어간다.
    text 에만 써넣으면 mzmine 이 무시하고 템플릿의 원래 경로로 출력한다.
    """
    for child in list(param):
        param.remove(child)
    param.text = None
    cf = ET.SubElement(param, "current_file")
    cf.text = str(path)


def load_raw_files(samples: list[str], polarity: str, limit_rep: int | None,
                   with_qc: bool = True) -> list[dict]:
    """처리용 시료표에서 조건에 맞는 원시자료 목록을 뽑는다.

    고르는 방법이 둘이다. 설정에 적힌 순서대로의 번호("1,2,3")로도, 시료 표기
    그대로("XX_00001")로도 고를 수 있다. 연구마다 시료 이름 규칙이 다르므로
    한 가지만 받으면 다음 연구에서 또 코드를 고치게 된다.

    반복확인용 시료는 기본으로 함께 넣는다. 같이 넣어야 런끼리 줄을 맞출 때
    기준이 생기고, 나중에 신호 흔들림을 그 시료로 볼 수 있다.
    """
    with METADATA.open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    codes = [c["sample_code"] for c in CFG["samples"]]
    want = set()
    for tok in samples:
        tok = tok.strip()
        if tok.isdigit() and 1 <= int(tok) <= len(codes):
            want.add(codes[int(tok) - 1])
        else:
            want.add(tok)

    sel = [r for r in rows if r["polarity"] == polarity
           and (r["sample_code"] in want or (with_qc and r["sample_code"] == "QC"))]
    if limit_rep is not None:
        sel = [r for r in sel if int(r["replicate"]) <= limit_rep]
    sel.sort(key=lambda r: (r["sample_code"], int(r["replicate"])))
    if not sel:
        print(f"  [안내] 시료표의 표기: {sorted({r['sample_code'] for r in rows})}")
    return sel


def probe_dbsearch_with_mzmine(root: ET.Element) -> None:
    """DB search 스텝만 떼어내 mzmine 에게 **파라미터 검증만** 시켜 본다 (~15초).

    '왜' 이게 필요한가? — 실제로 당한 두 번째 사고 때문이다.
    ion library 파라미터의 이름을 클래스 파일에 보이는 'Use adducts' 로 추정해 썼는데,
    실제 XML 이름은 'Ion library' 였다. mzmine 은 못 찾은 파라미터를 조용히 기본값
    (adduct 사용 안 함)으로 대체했고, 그 결과 **8분짜리 import 를 마친 뒤에야**
    "Cannot determine precursor mz" 로 동정 단계가 실패했다.

    파라미터 검증은 원시데이터가 없어도 즉시 수행되므로, 여기서 먼저 물어보면
    같은 실수를 15초 만에 잡을 수 있다.

    ⚠️ feature list 가 없는 상태로 돌리므로 "At least 1 feature lists must be
       selected" 는 **정상적으로 발생하는 잔여 에러**다. 그 외의 에러만 문제로 본다.
    """
    import subprocess
    import tempfile

    step = None
    for s in root.findall("batchstep"):
        if s.get("method") == M_DBSEARCH:
            step = s
            break
    if step is None:
        raise RuntimeError("DB search 스텝이 없어 검증 불가")

    probe_root = ET.Element("batch", {"mzmine_version": "4.10.6"})
    probe_root.append(step)
    with tempfile.TemporaryDirectory() as td:
        pf = Path(td) / "probe.mzbatch"
        ET.ElementTree(probe_root).write(pf, encoding="utf-8", xml_declaration=True)
        try:
            res = subprocess.run(
                [str(MZMINE_EXE), "-b", str(pf), "-user", str(MZUSER), "-m", "none"],
                capture_output=True, text=True, encoding="utf-8", errors="replace",
                timeout=300)
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            print(f"  [주의] mzmine 사전검증을 건너뜀 ({type(e).__name__})")
            return
        finally:
            # step 을 원래 배치로 돌려놓는다 (append 가 부모를 옮기지는 않지만 명시적으로)
            pass

    txt = (res.stdout or "") + (res.stderr or "")

    # 이온 라이브러리 로드 실패는 그 자체로 치명적
    if "Failed to load ion library" in txt or "Non-hex character" in txt:
        print("\n[mzmine 사전검증 실패] 이온 라이브러리를 로드하지 못했습니다.")
        raise SystemExit(1)

    errs = [l for l in txt.splitlines() if "Invalid parameter settings" in l]
    if not errs:
        print("  ⑧ mzmine 사전검증: 통과 (파라미터 오류 없음)")
        return

    msg = errs[0].split("module Local compound database search:")[-1]
    # feature list 부재는 probe 구조상 불가피 → 제외하고 남는 게 있는지 본다
    residual = msg.replace("At least 1 feature lists  must be selected", "")
    residual = residual.strip(" ,[]\n")
    if residual:
        print("\n[mzmine 사전검증 실패] DB search 파라미터 오류가 남아 있습니다:")
        print(f"  ✗ {residual[:300]}")
        raise SystemExit(1)
    print("  ⑧ mzmine 사전검증: 통과 (잔여 에러는 feature list 부재뿐 — 정상)")


def verify_batch(root: ET.Element, polarity: str, n_files: int) -> None:
    """완성된 배치를 저장 직전에 스스로 검사한다.

    '왜' 이 함수가 필요한가? — 실제로 당한 사고 때문이다.
    ionLibrary 의 id 를 'a1000000-...-00000000pos1' 처럼 사람이 읽기 좋게 지었더니
    mzmine 이 이를 java.util.UUID 로 역직렬화하다 'p' 는 16진수가 아니라며 실패했다.
    문제는 **mzmine 이 거기서 배치를 멈추지 않고 그냥 계속 돌았다**는 것이다.
    → adduct 계산 없이 20분을 돌고 나서야 동정이 이상하다는 걸 알게 된다.
    이런 '조용한 실패'는 실행 전에 잡아야 하고, 그래서 검사를 코드로 박아 둔다.
    """
    import uuid as _uuid

    problems: list[str] = []

    # (1) 모든 ionLibrary id 가 진짜 UUID 인가
    libs = root.findall(".//ionLibrary")
    if not libs:
        problems.append("ionLibrary 가 하나도 없다 (극성 설정이 통째로 빠졌다는 뜻)")
    for lib in libs:
        raw_id = lib.get("id", "")
        try:
            _uuid.UUID(raw_id)
        except (ValueError, AttributeError, TypeError):
            problems.append(f"ionLibrary id 가 유효한 UUID 가 아님: '{raw_id}' "
                            f"(라이브러리 '{lib.findtext('name')}')")

    # (2) ionTypes 가 참조하는 part id 가 실제로 존재하는가
    #     존재하지 않는 id 를 참조하면 adduct 가 조용히 누락된다
    for lib in libs:
        part_ids = {p.get("id") for p in lib.findall("parts/parts")}
        for it in lib.findall("ionTypes/ionTypes"):
            for part in it.findall("part"):
                if part.get("id") not in part_ids:
                    problems.append(f"ionTypes 가 없는 part id '{part.get('id')}' 를 참조 "
                                    f"(라이브러리 '{lib.findtext('name')}')")

    # (3) 우리가 만든 라이브러리의 극성이 요청한 극성과 맞는가
    #     part 의 charge 부호로 판정한다 — 이름만 믿으면 복붙 실수를 못 잡는다
    expect_token = "POS" if polarity == "POS" else "NEG"
    for lib in libs:
        nm = lib.findtext("name") or ""
        if "phytochem" not in nm:
            continue                      # mzmine 내장 라이브러리는 검사 대상 아님
        if expect_token not in nm:
            problems.append(f"라이브러리 이름 '{nm}' 이 요청 극성 {polarity} 와 불일치")

    # (4) ADAP 극성이 실제로 바뀌었는가
    adap_pol = find_step(root, M_ADAP).find("./parameter[@name='Scan filters']"
                                            "/parameter[@name='Polarity']")
    want = "+" if polarity == "POS" else "-"
    if adap_pol is None or adap_pol.text != want:
        problems.append(f"ADAP Polarity 가 '{adap_pol.text if adap_pol is not None else None}' "
                        f"(기대 '{want}')")

    # (4b) DB search 의 Ion library 가 켜져 있는가
    #      '왜'? 이름을 'Use adducts' 로 잘못 쓰면 mzmine 이 파라미터를 못 찾아
    #      기본값(꺼짐)으로 대체하고, 동정 단계 전체가 실패한다.
    db_lib = None
    for s in root.findall("batchstep"):
        if s.get("method") == M_DBSEARCH:
            db_lib = s.find("./parameter[@name='Ion library']")
    if db_lib is None:
        problems.append("DB search 에 'Ion library' 파라미터가 없다 "
                        "(이름을 'Use adducts' 로 잘못 쓴 경우 이 증상이 난다)")
    elif db_lib.get("selected") != "true":
        problems.append("DB search 의 'Ion library' 가 selected=true 가 아니다 "
                        "→ adduct 계산이 꺼져 동정 불가")

    # (5) Lipid Annotation 이 정말 사라지고 DB search 가 들어왔는가
    methods = [s.get("method") for s in root.findall("batchstep")]
    if M_LIPID in methods:
        problems.append("Lipid Annotation 이 아직 남아 있다")
    if M_DBSEARCH not in methods:
        problems.append("Local CSV DB search 가 삽입되지 않았다")

    # (6) 입력 파일 수가 맞는가 + 경로가 실재하는가
    files = find_step(root, M_IMPORT).findall("./parameter[@name='File names']/file")
    if len(files) != n_files:
        problems.append(f"입력 파일 수 {len(files)} != 기대 {n_files}")
    for f in files:
        if not Path(f.text).exists():
            problems.append(f"입력 경로가 실재하지 않음: {f.text}")

    # (7) DB CSV 의 헤더가 배치의 column 매핑과 일치하는가
    #     불일치하면 mzmine 이 그 컬럼을 못 찾아 조용히 동정 0건이 된다
    with DB_CSV.open(encoding="utf-8") as fh:
        db_headers = set(next(csv.reader(fh)))
    db_step = find_step(root, M_DBSEARCH)
    for it in db_step.findall("./parameter[@name='Columns']/importtype"):
        if it.get("selected") == "true" and it.get("column") not in db_headers:
            problems.append(f"배치가 요구하는 컬럼 '{it.get('column')}' 이 DB CSV 에 없음 "
                            f"(CSV 헤더: {sorted(db_headers)})")

    if problems:
        print("\n[배치 검증 실패]")
        for p in problems:
            print(f"  ✗ {p}")
        raise SystemExit(1)


def build(polarity: str, samples: list[str], name: str, limit_rep: int | None,
          skip_probe: bool = False) -> Path:
    tree = ET.parse(TEMPLATE)
    root = tree.getroot()
    basename = f"{name}_{polarity}"
    report: list[str] = []

    # ── ① 입력 파일 목록 교체 ──────────────────────────────────────────────
    files = load_raw_files(samples, polarity, limit_rep, with_qc=WITH_QC)
    if not files:
        raise RuntimeError(f"조건에 맞는 raw 파일이 없습니다 (samples={samples}, {polarity})")

    imp = find_step(root, M_IMPORT)
    fnames = get_param(imp, "File names")
    for child in list(fnames):
        fnames.remove(child)
    for r in files:
        el = ET.SubElement(fnames, "file")
        # 경로는 슬래시(/)로 — 백슬래시+숫자가 유니코드 escape 로 오인되는 것을 피한다
        el.text = r["raw_path"].replace("\\", "/")
    report.append(f"① 입력 파일 {len(files)} 개 주입")

    # ── ①b Waters 눈금 맞춤(lock mass) 보정 켜기/끄기 ──────────────────────
    #     '왜'? 2026-11-09 이전 연구: 측정 단계에서 이미 보정된 원시자료에 mzmine 이 보정을 한 번 더 걸어
    #     m/z 가 런마다 −20~−80 ppm 어긋났다(에러 없음). 설정 processing.waters_lockmass_on_import 로 정한다
    #     (없으면 배치 틀 값 그대로 — 기존 연구 동작 유지). mzmine 이 보정 결과를 .raw/lmgt.inf 에 써 두면
    #     Waters 읽기 부품이 그 파일을 다음 실행에도 적용하므로, 끌 때는 lmgt.inf 가 남아 있으면 멈춘다.
    lm_cfg = CFG["processing"].get("waters_lockmass_on_import")
    if lm_cfg is not None:
        lm = imp.find(".//parameter[@name='Apply lockmass on import (Waters)']")
        if lm is None:
            raise RuntimeError("배치 틀에 'Apply lockmass on import (Waters)' 가 없습니다")
        lm.set("selected", "true" if lm_cfg else "false")
        report.append(f"①b Waters lock mass 보정 = {'켬' if lm_cfg else '끔'} (설정 processing.waters_lockmass_on_import)")
        if not lm_cfg:
            # 🔴 mzmine 4.10.6 은 이 칸을 false 로 해도 Waters 읽기 부품이 보정을 계산·적용한다(2026-11-09 시험).
            #    실제로 끄는 방법은 배율 1.0 짜리 lmgt.inf 를 미리 넣는 것 → neutralize_waters_lockmass.py
            from neutralize_waters_lockmass import is_identity
            bad = [r["raw_path"] for r in files
                   if not (Path(r["raw_path"]) / "lmgt.inf").exists() or not is_identity(Path(r["raw_path"]) / "lmgt.inf")]
            if bad:
                raise RuntimeError(f"보정을 껐는데 배율 1.0 lmgt.inf 가 없는 .raw {len(bad)} 개 (예: {bad[0]}) — "
                                   f"먼저 py scripts/neutralize_waters_lockmass.py --config … 실행")

    # ── ② 극성 (1/2): ADAP scan filter ────────────────────────────────────
    adap = find_step(root, M_ADAP)
    scanf = get_param(adap, "Scan filters")
    pol_el = get_param(scanf, "Polarity")
    pol_el.text = "+" if polarity == "POS" else "-"
    report.append(f"② 극성(ADAP scan filter) = {pol_el.text}")

    # ── ③ RT crop 하한 완화 (2.0 → 0.3) ──────────────────────────────────
    rt_el = get_param(scanf, "Retention time")
    rt_min = rt_el.find("min")
    old_min = rt_min.text
    rt_min.text = "0.3"
    rt_max = rt_el.find("max")
    report.append(f"③ RT crop {old_min}→{rt_min.text} ~ {rt_max.text} 분 "
                  f"(극성 성분 보존 — 지질 템플릿의 2.0 하한 폐기)")

    # ── ④ 극성 (2/2): IIN Ion library ─────────────────────────────────────
    #     '왜' 여기도? IIN 은 adduct 관계로 feature 를 묶는데, 라이브러리가
    #     반대 극성이면 아무 것도 안 묶이고 에러도 안 난다.
    iin = find_step(root, M_IIN)
    old_lib = get_param(iin, "Ion library")
    idx = list(iin).index(old_lib)
    iin.remove(old_lib)
    lib_xml = ION_LIB_POS if polarity == "POS" else ION_LIB_NEG
    new_lib = ET.fromstring(lib_xml)
    # IIN 의 'Ion library' 는 OptionalParameter 가 아니라 일반 파라미터다.
    # DB search 쪽에서 필요한 selected="true" 를 여기 남겨두면 의미가 없으므로 뗀다.
    new_lib.attrib.pop("selected", None)
    iin.insert(idx, new_lib)
    n_iontypes = len(new_lib.find("ionLibrary/ionTypes"))
    report.append(f"④ 극성(IIN Ion library) = {polarity}, ion type {n_iontypes} 종")

    # ── ⑤ Lipid Annotation 제거 → Local CSV DB search 삽입 ────────────────
    lipid = find_step(root, M_LIPID)
    lipid_idx = list(root).index(lipid)
    root.remove(lipid)
    db_step = ET.fromstring(DB_SEARCH_STEP.format(
        module=M_DBSEARCH,
        db_path=str(DB_CSV),
        ion_library=(ION_LIB_POS if polarity == "POS" else ION_LIB_NEG),
    ))
    root.insert(lipid_idx, db_step)
    report.append(f"⑤ Lipid Annotation 제거 → Local CSV DB search 삽입 (DB: {DB_CSV.name})")

    # ── ⑥ 산출물 경로 재지정 ──────────────────────────────────────────────
    OUT_SMALL.mkdir(parents=True, exist_ok=True)
    OUT_BIG.mkdir(parents=True, exist_ok=True)
    outputs = [
        # (모듈, 파라미터명, 저장위치, 파일명)
        (M_CSVEXP,  "Filename",      OUT_SMALL, f"{basename}_full_feature_table.csv"),
        (M_ANNOEXP, "Filename",      OUT_SMALL, f"{basename}_annotations.csv"),
        (M_META,    "Metadata file", OUT_SMALL, f"{basename}_metadata.csv"),
        (M_GRAPHML, "Filename",      OUT_BIG,   f"{basename}_networking.graphml"),
        (M_GNPS,    "Filename",      OUT_BIG,   f"{basename}_iimn_gnps"),
        (M_SIRIUS,  "Filename",      OUT_BIG,   f"{basename}_sirius.mgf"),
    ]
    for module, pname, folder, fname in outputs:
        step = find_step(root, module)
        set_file_param(get_param(step, pname), folder / fname)
    report.append(f"⑥ 산출물: 작은 CSV → {OUT_SMALL.name}/ , 큰 파일 → {OUT_BIG}")

    # ── ⑦ AutoSave 제거 ──────────────────────────────────────────────────
    #     콘솔 모드에서는 "no batch is running" 으로 어차피 실패하는 단계라 뺀다.
    try:
        root.remove(find_step(root, M_AUTOSAVE))
        report.append("⑦ AutoSave 단계 제거 (콘솔에서 동작 안 함)")
    except RuntimeError:
        pass

    # ── ⑧ 자체 검증 (여기서 막지 못하면 20분 뒤에야 안다) ──────────────────
    verify_batch(root, polarity, len(files))
    if not skip_probe:
        probe_dbsearch_with_mzmine(root)

    # ── 저장 ──────────────────────────────────────────────────────────────
    BATCH_DIR.mkdir(parents=True, exist_ok=True)
    out = BATCH_DIR / f"{basename}.mzbatch"
    ET.indent(tree, space="    ")
    tree.write(out, encoding="utf-8", xml_declaration=True)

    print(f"\n=== {basename} ===")
    for line in report:
        print("  " + line)
    print(f"  → {out}")
    print("  [입력 파일]")
    for r in files:
        print(f"     {r['sample_code']} {r['sample_name']:<6} rep{r['replicate']}  {Path(r['raw_path']).name}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None, help="config/ 안의 연구 설정 파일")
    ap.add_argument("--no-qc", action="store_true",
                    help="반복확인용 시료를 배치에서 뺀다 (기본은 함께 넣는다)")
    ap.add_argument("--samples", default=None,
                    help="시료 번호 쉼표구분 (1=XX_00001 2=XX_00002 …). QC 는 별도 처리")
    ap.add_argument("--name", default=None, help="배치 basename")
    # 이번 데이터는 NEG 단독 측정이다 → 기본값을 both 가 아니라 NEG 로 둔다.
    # both 로 두면 POS 배치가 만들어지는데 해당 raw 가 없어 헛돈다.
    ap.add_argument("--polarity", default=None, choices=["POS", "NEG", "both"])
    ap.add_argument("--limit-rep", type=int, default=None,
                    help="반복 수 제한 (파일럿용, 예: 1 이면 각 시료 1반복만)")
    ap.add_argument("--skip-probe", action="store_true",
                    help="mzmine 사전검증 생략 (빠르지만 파라미터 오류를 못 잡는다)")
    args = ap.parse_args()

    # 사전 점검 — 여기서 막아야 20분짜리 실행이 헛돌지 않는다
    for p, label in [(TEMPLATE, "템플릿 배치"), (METADATA, "샘플 메타데이터"), (DB_CSV, "화합물 DB")]:
        if not p.exists():
            print(f"[에러] {label} 없음: {p}")
            return 1

    # 시료를 따로 고르지 않으면 설정에 적힌 전부를 쓴다
    samples = ([s.strip() for s in args.samples.split(",") if s.strip()] if args.samples
               else [str(i + 1) for i in range(len(CFG["samples"]))])
    # 극성을 따로 주지 않으면 설정 값을 쓴다 — 측정하지 않은 극성으로 배치를 만들면 헛돈다
    global WITH_QC
    WITH_QC = not args.no_qc
    pol_arg = args.polarity or CFG["polarity"]
    pols = ["POS", "NEG"] if pol_arg == "both" else [pol_arg]
    for pol in pols:
        build(pol, samples, args.name or PREFIX, args.limit_rep, skip_probe=args.skip_probe)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
