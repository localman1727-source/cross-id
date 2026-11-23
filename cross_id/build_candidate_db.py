# -*- coding: utf-8 -*-
"""
③ 후보 성분 목록 만들기 — 공개 천연물 자료에서 분류군 기준으로
================================================================
'왜' 시료에서 찾은 것으로 목록을 만들지 않는가?
  기기 제조사 소프트웨어가 바로 그 시료에서 찾아낸 목록으로 대조하면, 맞을 수밖에 없는
  것과 맞춰보게 되어 정확도가 낙관적으로 부풀고 목록의 근거를 논문에 쓸 수 없다.

'왜' 이 공개 자료(LOTUS)인가?
  천연물 자료 가운데 구조-생물종 짝마다 출처 논문이 붙어 있는 거의 유일한 공개 자료다.
  "이 성분을 왜 후보에 넣었나"가 논문 단위로 추적된다.
  출처: Rutz et al., eLife 2022;11:e70780 / 자료 Zenodo 10.5281/zenodo.19360665 (CC-BY-4.0)

'왜' 종이 아니라 속 단계인가?  ← 추측이 아니라 실측 결과다
  우리나라 특산종·희귀종은 공개 자료 수록이 0 인 경우가 흔하다. 종 단계로 목록을 만들면
  후보가 없어 조용히 빈 결과가 나온다. 속 단계는 "같은 속의 가까운 종에서 보고된 성분"
  이라는 가정 위에 서고, 그 가정은 taxon_rank 칸에 성분마다 기록돼 하류에서 구분된다.

'왜' 과(科) 단계는 기본값이 아닌가?
  과를 통째로 넣으면 같은 과의 큰 속(연구가 많이 된 작물 등)이 목록의 대부분을 차지해
  검색 공간만 넓어진다. 설정의 far_genera 에 그런 속을 적으면 근연 속 목록이 따로 나오고,
  그것을 실사용 후보로 삼는다.

⚠️ 이 목록으로 얻는 동정의 수준
  정확질량 대조이므로 잠정(분자식 수준)이다. 같은 분자식의 다른 구조는 원리적으로
  구별되지 않는다. 확정하려면 표준품을 같은 조건에서 돌려 머무는 시간까지 맞춰야 한다.

산출물 (runs/<study_id>/db/ 아래, 이름 앞에 연구 이름이 붙는다)
  *_candidates_genus.csv      속 단계
  *_candidates_family.csv     속+과 단계 (넓힌 판)
  *_candidates_relatives.csv  근연 속 목록  ★ 실사용
  *_ref_masses.csv            질량축 점검용 이론 질량 (극성에 맞춰 계산)
  *_candidates_check.json     검증 수치 전량 (규정 5-2, 반올림 없음)

실행: py scripts/build_candidate_db.py --config config/<연구>.json
"""
import sys
import io
import json
import re
import gzip
import urllib.request
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from study_config import load_study

CFG = load_study()
DB = CFG["db_dir"]
PREFIX = CFG["study_id"]
POLARITY = CFG["polarity"].upper()

# 원본 LOTUS 덤프 캐시.
# '왜' 홈 디렉토리인가? 90 MB 를 프로젝트 폴더(=보관소 동기화 대상)에 두면
# push 마다 따라다닌다. 홈 아래면 계정명이 경로에 박히지 않으면서(규정 3-1)
# 동기화 대상도 아니다. 다른 프로젝트에서도 재사용된다.
CACHE = Path.home() / ".npdb_cache"
LOTUS_URL = ("https://zenodo.org/api/records/19360665/files/"
             "260413_frozen_metadata.csv.gz/content")
LOTUS_GZ = CACHE / "lotus_260413_frozen_metadata.csv.gz"
LOTUS_CITATION = {
    "dataset": "The LOTUS Initiative for Open Natural Products Research: "
               "frozen dataset union wikidata (with metadata)",
    "zenodo_doi": "10.5281/zenodo.19360665",
    "frozen_date": "2026-04-13",
    "license": "CC-BY-4.0",
    "paper_doi": "10.7554/eLife.70780",
    "downloaded_sha256": "(빌드 시 계산)",
}

# ── 대상 분류군 (설정 파일에서 읽는다) ─────────────────────────────────────
#   국명은 시료 라벨 표기가 아니라 국가표준식물목록 표준 국명을 쓴다. 라벨이 표준명과
#   다른 경우가 흔하고(오기·옛 이름), 그대로 두면 검색이 빗나간다.
#
# 🔴 과(科) 이름이 둘로 돌아다니는 종이 많다 — 한쪽만 찾으면 조용히 절반을 놓친다
#   분류가 개정되면 옛 문헌과 일부 공개 자료는 여전히 옛 과명으로 싣는다. 공개 자료는
#   원논문 표기를 따르므로 두 이름을 모두 찾아야 한다. 설정의 family_aliases 가 그 자리다.
TAXA = [(s_.get("sample_code", ""), s_.get("collect_id", ""), s_.get("name_kr", ""),
         s_.get("scientific_name", ""), s_.get("genus", ""), s_.get("family", ""),
         s_.get("note", "")) for s_ in CFG["samples"]]

FAMILY_ALIASES = CFG.get("family_aliases") or []

GENERA = sorted({t[4] for t in TAXA if t[4]})
FAMILIES = sorted({f for f in (set(FAMILY_ALIASES) | {t[5] for t in TAXA}) if f})

# ── 단일동위원소질량 (IUPAC/CIAAW) ──────────────────────────────────────────
#   '왜' LOTUS 의 structure_exact_mass 를 그대로 안 쓰나?
#   ① 남의 계산을 검증 없이 신뢰하면 그 오류가 조용히 전파된다.
#   ② 우리 질량축 검증(diag_mass_offset·crosscheck_vendor_mz)이 이 값을 기준으로
#      삼는데, 외부 값을 쓰면 "외부 축으로 외부 축을 검증"하는 순환이 된다.
#   대신 계산값을 LOTUS 값과 **대조**해 파서 자체를 검증한다(둘이 맞으면 양쪽 다 옳다).
ATOMIC = {
    "C": 12.0, "H": 1.0078250319, "N": 14.0030740052, "O": 15.9949146221,
    "S": 31.97207069, "P": 30.97376151,
    "F": 18.99840320, "Cl": 34.96885271, "Br": 78.9183376, "I": 126.904473,
    "Na": 22.98976928, "K": 38.96370649, "Mg": 23.98504170, "Ca": 39.96259098,
    "Fe": 55.93494210, "Zn": 63.92914220, "Mn": 54.93804510, "Cu": 62.92959750,
    "Si": 27.97692653, "Se": 79.9165218, "B": 11.0093054, "Al": 26.98153853,
}
ELECTRON = 0.00054857990927
PROTON = 1.00727645209073        # H 원자량 − 전자질량

FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def exact_mass(formula):
    """분자식 → 단일동위원소 중성질량. 해석 실패 시 None.

    정규식이 소비한 글자 수를 원본 길이와 대조한다.
    '왜'? 괄호·전하표기·소문자 시작 같은 이상 문자열이 있으면 findall 은
    조용히 일부만 먹고 넘어가 **틀린 질량이 에러 없이** 만들어진다.
    """
    f = str(formula).strip()
    if not f or f.lower() == "nan":
        return None
    consumed, total = 0, 0.0
    for m in FORMULA_RE.finditer(f):
        sym, cnt = m.group(1), m.group(2)
        if sym not in ATOMIC:
            return None
        total += ATOMIC[sym] * (int(cnt) if cnt else 1)
        consumed += len(m.group(0))
    if consumed != len(f):
        return None
    return total


def count_element(formula, elem):
    """분자식 안 특정 원소 개수 (카로티노이드 carotene/xanthophyll 판별용)."""
    for m in FORMULA_RE.finditer(str(formula)):
        if m.group(1) == elem:
            return int(m.group(2)) if m.group(2) else 1
    return 0


# ── NPClassifier class → refine_annotations.RT_WINDOW 키 매핑 ───────────────
#   🔴 설계 원칙: **화학적으로 명확한 것만 매핑한다.**
#
#   '왜' 전부 매핑하지 않나?
#     refine_annotations.py 는 RT_WINDOW 에 없는 클래스를 (0.0, 99.0) 으로 처리한다
#     = RT 증거가 무력화될 뿐 화합물이 버려지지는 않는다(IIN·공용출 증거로 Tier1 가능).
#     반대로 **틀린 RT 창을 주면 실존 화합물이 조용히 Tier3 로 강등된다.**
#     즉 "모르면 비워 두는 것"이 안전한 방향이고, 추측 창을 새로 만드는 것이 위험하다.
#     (구 DB 는 전 항목이 미매핑이라 RT 증거가 통째로 꺼져 있었다 — 이것도 문제지만,
#      해법은 '추측으로 채우기'가 아니라 '근거 있는 것만 채우기'다.)
#
#   미매핑 화합물은 `Unmapped_<pathway>` 로 남겨 개수를 추적한다.
NPC_CLASS_MAP = {
    # 플라보노이드 — 하위 골격이 RT 를 잘 가른다
    "Flavonols": "Flavonoid_flavonol",
    "Flavones": "Flavonoid_flavone",
    "Flavanones": "Flavonoid_flavanone",
    "Dihydroflavonols": "Flavonoid_flavanonol",
    "Flavan-3-ols": "Flavonoid_flavan3ol",
    "Flavanols": "Flavonoid_flavan3ol",
    "Isoflavones": "Flavonoid_isoflavone",
    "Isoflavonoids": "Flavonoid_isoflavone",
    "Proanthocyanins": "Flavonoid_proanthocyanidin",
    "Proanthocyanidins": "Flavonoid_proanthocyanidin",
    # 페놀산 — 극성, 이른 용출
    "Cinnamic acids and derivatives": "Phenolic_acid",
    "Simple phenolic acids": "Phenolic_acid",
    "Benzoic acids and derivatives": "Phenolic_acid",
    "Hydroxybenzoic acids": "Phenolic_acid",
    "Hydroxycinnamic acids": "Phenolic_acid",
    # 스틸벤 · 쿠마린 · 리그난
    "Monomeric stilbenes": "Stilbene",
    "Stilbenes": "Stilbene",
    "Oligomeric stilbenes": "Stilbene",
    "Simple coumarins": "Coumarin",
    "Coumarins": "Coumarin",
    "Furanocoumarins": "Coumarin",
    "Pyranocoumarins": "Coumarin",
    "Arylnaphthalene and aryltetralin lignans": "Lignan",
    "Dibenzylbutane lignans": "Lignan",
    "Dibenzylbutyrolactone lignans": "Lignan",
    "Furanoid lignans": "Lignan",
    "Furofuranoid lignans": "Lignan",
    "Lignans": "Lignan",
    # 트리테르페노이드 — 무극성, 늦은 용출
    "Dammarane and Protostane triterpenoids": "Triterpenoid",
    "Lupane triterpenoids": "Triterpenoid",
    "Oleanane triterpenoids": "Triterpenoid",
    "Ursane and Taraxastane triterpenoids": "Triterpenoid",
    "Taraxerane triterpenoids": "Triterpenoid",
    "Lanostane, Tirucallane and Euphane triterpenoids": "Triterpenoid",
    "Cycloartane triterpenoids": "Triterpenoid",
    "Hopane triterpenoids": "Triterpenoid",
    "Friedelane triterpenoids": "Triterpenoid",
    "Cucurbitane triterpenoids": "Triterpenoid",
    # 스테롤
    "Stigmastane steroids": "Sterol",
    "Cholestane steroids": "Sterol",
    "Ergostane steroids": "Sterol",
    "Steroids": "Sterol",
    # 1차 대사체 — 극성, void 근처
    "Monosaccharides": "Primary_metabolite",
    "Disaccharides": "Primary_metabolite",
    "Oligosaccharides": "Primary_metabolite",
    "Aminoacids": "Primary_metabolite",
    "Amino acids": "Primary_metabolite",
    "Sugar alcohols": "Primary_metabolite",
    "Tricarboxylic acids": "Primary_metabolite",
    # 카로티노이드는 분자식의 산소 유무로 하위분류 (아래 assign_class 에서 처리)
}

# 카로티노이드는 class 문자열이 "Carotenoids (C40, β-β)" 처럼 변형이 많아 접두 매칭
NPC_CAROTENOID_PREFIX = "Carotenoid"


def assign_class(npc_class, npc_pathway, formula):
    """NPClassifier 분류 → refine_annotations 의 molecular_class.

    매핑이 명확하지 않으면 `Unmapped_<pathway>` 를 돌려준다
    (RT_WINDOW 에 없으므로 RT 증거가 비활성 = 보수적 동작).
    """
    c = str(npc_class).strip()
    p = str(npc_pathway).strip().replace(" ", "_").replace("$", "and")

    if c and c.lower() != "nan":
        if c.startswith(NPC_CAROTENOID_PREFIX):
            # C40 탄화수소 = carotene / 산소 함유 = xanthophyll.
            # '왜' 이렇게 가르나? 이 구분이 곧 극성 차이이고 RT 를 가르는 실제 요인이다.
            return ("Carotenoid_xanthophyll" if count_element(formula, "O") > 0
                    else "Carotenoid_carotene")
        if c in NPC_CLASS_MAP:
            return NPC_CLASS_MAP[c]

    if not p or p.lower() == "nan":
        return "Unmapped_unknown"
    return f"Unmapped_{p}"


# ── 이온화 타당성 필터 (NEG 전용) — 규정 ⑩ 의 지질 타당성 필터에 대응 ──────
#   '왜' tier/decoy 위에 이 관문이 따로 필요한가?
#     tier 는 "질량이 맞고 증거가 있는가"만 본다. **그 화합물이 이 극성에서
#     애초에 이온이 되는가**는 묻지 않는다. 그래서 NEG 단독 측정인데도
#     카로티노이드처럼 [M-H]- 를 만들 수 없는 물질이 후보로 남는다.
#     (지질 쪽에서 tier 가 'PC O-19:0' 같은 유령을 Tier1 으로 승격시킨 것과 같은 구조의 문제.)
#
#   2026-08-06 실측이 이 필터의 근거다:
#     LOTUS 속 DB 902종 중 카로티노이드 49종이 들어 있었고,
#     본 시료 종에서 직접 보고된 122종 중 48종(39%)이 카로티노이드였다.
#     이들은 NEG 에서 Tier1 도달 0종 — 즉 검색 공간만 넓히고 기여는 없었다.
#
#   🔴 제거 vs 플래그 구분 (이전 작업의 원칙 계승)
#     · 제거: 물리적으로 불가능한 것만. 산소·질소가 하나도 없으면
#             내놓을 양성자가 없어 [M-H]- 자체가 성립하지 않는다.
#     · 플래그: 경험적으로 약한 것(카로티노이드 등). 지우지 않고 표시만 해
#             하류에서 선택할 수 있게 둔다 — 근거가 경험적이면 되돌릴 수 있어야 한다.
NEG_WEAK_CLASSES = ("Carotenoid_carotene", "Carotenoid_xanthophyll")


def neg_ionization_flags(formula, molecular_class):
    """NEG 모드 이온화 타당성. 반환 (removable, weak_flag, reason)."""
    n_o = count_element(formula, "O")
    n_n = count_element(formula, "N")
    n_s = count_element(formula, "S")
    n_p = count_element(formula, "P")

    # 헤테로원자가 전무 = 산성 양성자 공여 불가 → [M-H]- 물리적으로 불가능
    if (n_o + n_n + n_s + n_p) == 0:
        return True, True, "no_heteroatom_cannot_deprotonate"

    if molecular_class in NEG_WEAK_CLASSES:
        return False, True, "carotenoid_poor_neg_ionization"

    return False, False, ""


def load_lotus():
    """LOTUS 덤프를 캐시에서 읽고, 없으면 내려받는다."""
    CACHE.mkdir(parents=True, exist_ok=True)
    if not LOTUS_GZ.exists():
        print(f"[다운로드] {LOTUS_URL}")
        print("  (90 MB · 최초 1회만. 이후 캐시 재사용)")
        urllib.request.urlretrieve(LOTUS_URL, LOTUS_GZ)
    import hashlib
    h = hashlib.sha256()
    with open(LOTUS_GZ, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    LOTUS_CITATION["downloaded_sha256"] = h.hexdigest()

    cols = ["structure_inchikey", "structure_molecular_formula", "structure_exact_mass",
            "structure_nameTraditional", "structure_smiles_2D", "structure_xlogp",
            "structure_taxonomy_npclassifier_01pathway",
            "structure_taxonomy_npclassifier_03class",
            "organism_taxonomy_06family", "organism_taxonomy_08genus",
            "organism_taxonomy_09species", "reference_doi"]
    return pd.read_csv(LOTUS_GZ, usecols=cols, low_memory=False)


def build_db(sub, rank_label, check):
    """분류군 부분집합 → mzmine 입력 DB (구조 단위로 축약 + provenance).

    rank_label: 'genus' | 'genus+family' — 어느 범위로 모았는지 기록용
    """
    # 🔴 인덱스를 반드시 초기화한다.
    #    '왜'? 속+과 판은 pd.concat 으로 만들어져 **인덱스 라벨이 중복**된다.
    #    그 상태에서 Series.where(cond, other_series) 를 쓰면 라벨 정렬이 어긋나
    #    엉뚱한 행에 값이 들어간다 — 에러 없이 조용히 틀리는 유형이다
    #    (2026-08-06 실제로 걸림: 이름 결측 135종이 NaN 인 채로 통과해
    #     하류 sorted() 가 float 와 str 을 비교하며 죽었다).
    sub = sub.reset_index(drop=True).copy()

    # 화합물 이름: 관용명 우선, 없으면 InChIKey (이름 없는 구조도 질량 매칭은 가능하다).
    # 인덱스 정렬을 아예 타지 않도록 numpy 배열로 처리한다.
    nm = sub["structure_nameTraditional"].astype(str).str.strip().to_numpy()
    ik = sub["structure_inchikey"].astype(str).to_numpy()
    sub["cmpd_name"] = [
        k if (not str(n) or str(n).lower() == "nan") else n for n, k in zip(nm, ik)
    ]

    # 구조(InChIKey) 단위로 축약. 같은 구조가 여러 종·문헌에서 보고되므로
    # 그 근거들을 모아 provenance 로 남긴다.
    def agg_taxa(g):
        sp = sorted({s for s in g.dropna().astype(str) if s and s != "nan"})
        return "|".join(sp[:8]) + ("|…" if len(sp) > 8 else "")

    grp = sub.groupby("structure_inchikey", sort=False)
    out = pd.DataFrame({
        "inchikey": grp.size().index,
        "name": grp["cmpd_name"].first().to_numpy(),
        "formula": grp["structure_molecular_formula"].first().to_numpy(),
        "lotus_exact_mass": grp["structure_exact_mass"].first().to_numpy(),
        "xlogp": grp["structure_xlogp"].first().to_numpy(),
        "npc_pathway": grp["structure_taxonomy_npclassifier_01pathway"].first().to_numpy(),
        "npc_class": grp["structure_taxonomy_npclassifier_03class"].first().to_numpy(),
        "source_species": grp["organism_taxonomy_09species"].apply(agg_taxa).to_numpy(),
        "source_genera": grp["organism_taxonomy_08genus"].apply(agg_taxa).to_numpy(),
        "source_families": grp["organism_taxonomy_06family"].apply(agg_taxa).to_numpy(),
        "n_references": grp["reference_doi"].nunique().to_numpy(),
        "reference_dois": grp["reference_doi"].apply(agg_taxa).to_numpy(),
    })

    # taxon_rank — 이 화합물이 '본 시료 종'에서 직접 보고된 것인지,
    #              같은 속/과의 근연종에서 온 것인지. 하류 해석에서 반드시 구분해야 한다.
    # 학명에서 명명자를 떼고 두 낱말(속 + 종소명)로만 맞춘다.
    # '왜'? 설정에는 "Genus species (Author) Author" 처럼 명명자까지 적는 것이 맞지만,
    # 공개 자료는 명명자를 빼고 싣는다. 그대로 대조하면 한 건도 안 맞으면서 에러도 안 난다 —
    # 종 단계 근거가 조용히 0 이 되고 속 단계 가정만 남는다.
    def _binomial(name):
        parts = str(name).replace("×", " ").split()
        return " ".join(parts[:2]) if len(parts) >= 2 else str(name)

    target_species = {_binomial(t[3]) for t in TAXA if t[3]}
    sp_hit = (sub[sub["organism_taxonomy_09species"].map(_binomial).isin(target_species)]
              .groupby("structure_inchikey").size())
    ge_hit = (sub[sub["organism_taxonomy_08genus"].isin(GENERA)]
              .groupby("structure_inchikey").size())
    out["taxon_rank"] = "family"
    out.loc[out["inchikey"].isin(ge_hit.index), "taxon_rank"] = "genus"
    out.loc[out["inchikey"].isin(sp_hit.index), "taxon_rank"] = "species"

    # ── 질량 재계산 + LOTUS 값과 대조 (파서 검증) ──────────────────────────
    out["neutral_mass"] = out["formula"].map(exact_mass)
    n_before = len(out)
    failed = out[out["neutral_mass"].isna()]
    out = out.dropna(subset=["neutral_mass"]).copy()

    out["diff_vs_lotus_mDa"] = (out["neutral_mass"] - out["lotus_exact_mass"]) * 1000
    bad = out[out["diff_vs_lotus_mDa"].abs() > 1.0]

    out["molecular_class"] = [
        assign_class(c, p, f) for c, p, f in
        zip(out["npc_class"], out["npc_pathway"], out["formula"])
    ]

    # ── 검증 수치 (규정 5-2: 반올림 없이 raw) ──────────────────────────────
    RT_MAPPED = ~out["molecular_class"].str.startswith("Unmapped_")
    key = f"db_{rank_label.replace('+', '_')}"
    check[key] = {
        "n_structures_before_parse": int(n_before),
        "n_parse_failed": int(len(failed)),
        "parse_failed_formulas": failed["formula"].astype(str).tolist()[:20],
        "n_structures": int(len(out)),
        "n_unique_formulas": int(out["formula"].nunique()),
        "n_unique_masses_1mDa": int(out["neutral_mass"].round(3).nunique()),
        "diff_vs_lotus_mDa_median": float(out["diff_vs_lotus_mDa"].median()),
        "diff_vs_lotus_mDa_max_abs": float(out["diff_vs_lotus_mDa"].abs().max()),
        "n_diff_over_1mDa": int(len(bad)),
        "taxon_rank_counts": {k: int(v) for k, v in
                              out["taxon_rank"].value_counts().items()},
        "n_rt_window_mapped": int(RT_MAPPED.sum()),
        "n_rt_window_unmapped": int((~RT_MAPPED).sum()),
        "molecular_class_counts": {k: int(v) for k, v in
                                   out["molecular_class"].value_counts().items()},
    }

    # ── NEG 이온화 타당성 플래그 ───────────────────────────────────────────
    flags = [neg_ionization_flags(f, c)
             for f, c in zip(out["formula"], out["molecular_class"])]
    out["neg_removable"] = [x[0] for x in flags]
    out["neg_weak"] = [x[1] for x in flags]
    out["neg_flag_reason"] = [x[2] for x in flags]
    check[key]["neg_ionization"] = {
        "n_removable_no_heteroatom": int(out["neg_removable"].sum()),
        "n_weak_total": int(out["neg_weak"].sum()),
        "reason_counts": {k: int(v) for k, v in
                          out.loc[out["neg_weak"], "neg_flag_reason"].value_counts().items()},
    }

    # ── 최종 가드 — 이름 결측은 하류에서 TypeError 로 터진다 ─────────────────
    #   '왜' assert 인가? 조용히 넘어가면 refine_annotations 의 sorted() 가
    #   float(NaN) 와 str 을 비교하다 죽는다. 여기서 막는 편이 진단이 쉽다.
    n_null_name = int(out["name"].isna().sum()) + int((out["name"].astype(str) == "").sum())
    if n_null_name:
        raise ValueError(f"[중단] 이름이 빈 화합물 {n_null_name} 종 — cmpd_name 대체 로직 확인")
    check[key]["n_null_name_after_guard"] = 0

    cols = ["name", "formula", "neutral_mass", "molecular_class",
            "inchikey", "taxon_rank", "source_species", "source_genera",
            "source_families", "n_references", "reference_dois",
            "npc_pathway", "npc_class", "xlogp",
            "neg_removable", "neg_weak", "neg_flag_reason",
            "lotus_exact_mass", "diff_vs_lotus_mDa"]
    return out[cols].sort_values("name").reset_index(drop=True), bad


def main():
    check = {"lotus_source": LOTUS_CITATION,
             "target_taxa": [{"sample": t[0], "id": t[1], "korean_name": t[2],
                              "scientific_name": t[3], "genus": t[4],
                              "family": t[5], "note": t[6]} for t in TAXA]}

    df = load_lotus()
    print(f"LOTUS 전체 {len(df):,} 행 (구조-생물종-문헌 삼중항)")
    print(f"  대상 속 {GENERA}")
    print(f"  대상 과 {FAMILIES}\n")

    gen_sub = df[df["organism_taxonomy_08genus"].isin(GENERA)]
    fam_sub = df[df["organism_taxonomy_06family"].isin(FAMILIES)]

    # ── 종 단계 커버리지 실측 (설계 근거를 매 빌드마다 재확인) ──────────────
    print("=== 분류 단계별 커버리지 (고유 구조 수) ===")
    print(f"{'시료':<4}{'국명':<14}{'학명':<26}{'종':>6}{'속':>7}{'과':>7}")
    print("-" * 66)
    cov = []
    for code, sid, kr, sci, genus, fam, note in TAXA:
        n_sp = df[df["organism_taxonomy_09species"] == sci]["structure_inchikey"].nunique()
        n_ge = df[df["organism_taxonomy_08genus"] == genus]["structure_inchikey"].nunique()
        n_fa = df[df["organism_taxonomy_06family"] == fam]["structure_inchikey"].nunique()
        print(f"{code:<4}{kr:<14}{sci:<26}{n_sp:>6,}{n_ge:>7,}{n_fa:>7,}")
        cov.append({"sample": code, "scientific_name": sci,
                    "n_species": int(n_sp), "n_genus": int(n_ge), "n_family": int(n_fa)})
    check["taxon_coverage"] = cov

    # ── DB 두 판 빌드 ──────────────────────────────────────────────────────
    print("\n=== DB 빌드 ===")
    db_gen, bad_gen = build_db(gen_sub, "genus", check)
    db_fam, bad_fam = build_db(pd.concat([gen_sub, fam_sub]).drop_duplicates(),
                               "genus+family", check)

    for label, d, bad in [("속 단계", db_gen, bad_gen),
                          ("속+과 단계", db_fam, bad_fam)]:
        k = "db_genus" if label == "속 단계" else "db_genus_family"
        c = check[k]
        print(f"\n[{label}]")
        print(f"  구조 {c['n_structures']:,} 종 "
              f"(고유 분자식 {c['n_unique_formulas']:,} · 고유질량 {c['n_unique_masses_1mDa']:,})")
        print(f"  ⚠️ MS1 매칭의 실효 크기는 '구조 수'가 아니라 '고유 질량 수'다 "
              f"— 이성질체는 질량이 같아 구별되지 않는다")
        print(f"  분자식 파싱 실패 {c['n_parse_failed']} 종")
        print(f"  질량 계산 vs LOTUS: 중앙값 {c['diff_vs_lotus_mDa_median']:.6f} mDa · "
              f"최대 {c['diff_vs_lotus_mDa_max_abs']:.6f} mDa · 1 mDa 초과 {c['n_diff_over_1mDa']} 종")
        if len(bad):
            print("  ⚠️ 1 mDa 넘게 어긋난 상위 5:")
            print(bad.nlargest(5, "diff_vs_lotus_mDa")[
                ["name", "formula", "neutral_mass", "lotus_exact_mass",
                 "diff_vs_lotus_mDa"]].to_string(index=False))
        print(f"  분류군 출처: {c['taxon_rank_counts']}")
        print(f"  RT 창 매핑됨 {c['n_rt_window_mapped']:,} / "
              f"미매핑 {c['n_rt_window_unmapped']:,} "
              f"({c['n_rt_window_unmapped']/c['n_structures']*100:.1f}%는 RT 증거 비활성)")
        ni = c["neg_ionization"]
        print(f"  NEG 이온화: 제거대상(헤테로원자 0) {ni['n_removable_no_heteroatom']:,} · "
              f"약함 플래그 총 {ni['n_weak_total']:,} {ni['reason_counts']}")

    # ── 저장 ───────────────────────────────────────────────────────────────
    #   🔴 BOM 없이 저장. utf-8-sig 로 쓰면 첫 컬럼명이 '﻿name' 이 되어
    #      mzmine 이 'name' 컬럼을 못 찾는다 (2026-07-31 실제로 걸린 함정).
    db_gen.to_csv(DB / f"{PREFIX}_candidates_genus.csv", index=False, encoding="utf-8")
    db_fam.to_csv(DB / f"{PREFIX}_candidates_family.csv", index=False, encoding="utf-8")

    # NEG 이온화 타당성 적용판 — 이번 측정(NEG 단독)의 실사용 후보
    #   '왜' 별도 파일로 두나? 원본을 지우면 "무엇이 왜 빠졌는지" 추적이 끊긴다.
    #   두 파일을 나란히 두면 필터의 효과를 언제든 되돌려 검증할 수 있다(규정 11-4 정신).
    # 이 걸러내기는 음이온 측정에서만 뜻이 있다. 양이온에서는 헤테로원자가 없어도
    # 라디칼 양이온으로 뜨므로 제외 근거가 성립하지 않는다 — 파일은 남기되 이름으로 구분한다.
    db_gen_neg = db_gen[~db_gen["neg_removable"] & ~db_gen["neg_weak"]].copy()
    db_gen_neg.to_csv(DB / f"{PREFIX}_candidates_genus_negfilter.csv", index=False, encoding="utf-8")
    if POLARITY != "NEG":
        print("  ⚠ 위 걸러내기 결과는 이번 측정(양이온)에 쓰지 않는다 — 양이온에서는 "
              "헤테로원자가 없어도 이온이 되므로 제외 근거가 성립하지 않는다")
    check["db_genus_negfilt"] = {
        "n_structures": int(len(db_gen_neg)),
        "n_unique_masses_1mDa": int(db_gen_neg["neutral_mass"].round(3).nunique()),
        "n_dropped": int(len(db_gen) - len(db_gen_neg)),
        "dropped_reason_counts": {k: int(v) for k, v in
                                  db_gen.loc[db_gen["neg_weak"], "neg_flag_reason"]
                                  .value_counts().items()},
        "taxon_rank_counts": {k: int(v) for k, v in
                              db_gen_neg["taxon_rank"].value_counts().items()},
    }
    cn = check["db_genus_negfilt"]
    print(f"\n[속 단계 · NEG 이온화 필터 적용]")
    print(f"  {len(db_gen):,} → {cn['n_structures']:,} 종 "
          f"(제외 {cn['n_dropped']:,}) · 고유질량 {cn['n_unique_masses_1mDa']:,}")
    print(f"  제외 사유: {cn['dropped_reason_counts']}")
    print(f"  분류군 출처: {cn['taxon_rank_counts']}")

    # 질량축 검증용 양이온 기준 (속 단계 기준, 분자식 중복 제거)
    # 🔴 이전 원본은 NEG 측정이라 [M-H]- 를 썼다. 이번 측정은 POS 단독이므로
    #    같은 식을 그대로 두면 파일 이름만 pos 이고 계산은 neg 인 상태가 된다 —
    #    질량축 점검이 통째로 어긋나면서도 에러는 나지 않는다.
    NA_ION = 22.98922070009073      # Na 원자량 − 전자질량
    NH4_ION = 18.03382554809073     # NH4 − 전자질량
    HCOO_ION = 44.99820285790927    # 포름산 이온 − 전자질량 규약 반영
    base = db_gen["neutral_mass"]
    if POLARITY == "POS":
        cols = {"mz_M+H": base + PROTON, "mz_M+Na": base + NA_ION, "mz_M+NH4": base + NH4_ION}
        sort_by = "mz_M+H"
    else:
        # 역상 이동상에 포름산을 쓰는 것이 표준이라, 중성 배당체·당류는 [M-H]- 가 약하고
        # 포름산이 붙은 형태로 주로 뜬다. 빼면 배당체를 놓친다.
        cols = {"mz_M-H": base - PROTON, "mz_M+HCOO": base + HCOO_ION}
        sort_by = "mz_M-H"
    ref = (pd.DataFrame({"formula": db_gen["formula"], "neutral_mass": base, **cols})
           .drop_duplicates("formula").sort_values(sort_by))
    ref.to_csv(DB / f"{PREFIX}_ref_masses.csv", index=False, encoding="utf-8")
    check["n_ref_masses"] = int(len(ref))

    # ── 근연 속 우선 목록 (이 프로젝트의 실사용 후보) ──────────────────────
    # 🔴 '왜' 과 단계를 그대로 쓰지 않는가?
    #   과 단계 1,068종 가운데 Cannabis 524 · Humulus 497 = 약 95%다. 두 속은 같은 과지만
    #   칸나비노이드·쓴맛 수지처럼 이 나무에 있을 이유가 없는 성분을 대량으로 들고 온다.
    #   그 목록으로 맞추면 후보는 늘고 맞을 확률은 떨어진다 — 검색 공간만 넓어진다.
    # '왜' 속 단계만으로는 안 되는가?
    #   대상 속 는 LOTUS 수록이 5종뿐이다. 그 목록으로는 거의 아무것도 못 맞춘다.
    #   즉 이 종은 "문헌에 기록이 없는 쪽"이지 "성분이 없는 쪽"이 아니다.
    # → 타협: 같은 과 안에서도 **계통적으로 가까운 속**만 남긴다.
    #   가까운 속(genus)과 옛 분류에서 같은 과였던 속을 함께 쓴다.
    FAR_GENERA = set(CFG.get("far_genera") or ())   # 같은 과지만 화학이 동떨어진 속

    def _has_near_genus(v):
        gs = {x.strip() for x in str(v).split("|")} if pd.notna(v) else set()
        return bool(gs - FAR_GENERA)

    db_rel = db_fam[db_fam["source_genera"].apply(_has_near_genus)].copy()
    db_rel.to_csv(DB / f"{PREFIX}_candidates_relatives.csv", index=False, encoding="utf-8")
    check["db_relatives"] = {
        "n_structures": int(len(db_rel)),
        "n_unique_formula": int(db_rel["formula"].nunique()),
        "excluded_genera": sorted(FAR_GENERA),
        "rationale": "같은 과이지만 화학적으로 동떨어진 속(Cannabis·Humulus)을 제외한 근연 속 목록",
    }
    print("")
    print("[근연 속 목록 — 이 프로젝트 실사용]")
    print(f"  {len(db_fam):,} → {len(db_rel):,} 종 (제외 속: {', '.join(sorted(FAR_GENERA))})")
    print(f"  고유 분자식 {db_rel['formula'].nunique():,}")

    with open(DB / f"{PREFIX}_candidates_check.json", "w", encoding="utf-8") as f:
        json.dump(check, f, ensure_ascii=False, indent=2)

    print(f"\n저장:")
    print(f"  {DB / (PREFIX + '_candidates_genus.csv')}      ({len(db_gen):,} 종)")
    print(f"  {DB / (PREFIX + '_candidates_family.csv')}     ({len(db_fam):,} 종)")
    print(f"  {DB / (PREFIX + '_candidates_relatives.csv')}  ★ 실사용")
    print(f"  {DB / (PREFIX + '_ref_masses.csv')}            ({len(ref):,} 질량 · {POLARITY})")
    print(f"  {DB / (PREFIX + '_candidates_check.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
