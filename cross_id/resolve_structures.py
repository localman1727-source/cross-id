# -*- coding: utf-8 -*-
"""
이름 → 구조 지문(InChIKey) 붙이기 — 두 도구가 같은 물질을 다른 이름으로 부를 때 가려내기
=========================================================================================
'왜' 이름이 아니라 구조 지문으로 비교하나?
  같은 물질도 사전마다 이름이 다르다(Scopoletin = 6-Hydroxy-7-methoxycoumarin,
  L-Tryptophan = (2S)-2-ammonio-3-(1H-indol-3-yl)propanoate). 이름 글자를 비교하면
  이런 동의어가 모두 "충돌"로 잡힌다. InChIKey 는 구조에서 계산되는 지문이라 이름과 무관하다.
    앞 14자  = 원자 연결(뼈대)
    가운데  = 입체·동위원소
  앞 14자가 같으면 같은 뼈대(입체만 다를 수 있음), 전체가 같으면 같은 구조다.

'왜' 지문을 먼저 파일에서 찾고 PubChem 은 나중인가?
  두 도구가 이름을 가져온 원본에 이미 지문이 있다 — MS-DIAL 은 공개 사전(.msp)의 INCHIKEY 줄,
  mzmine 은 LOTUS 에서 만든 후보 목록(candidates_*.csv)의 inchikey 칸. 원본에서 찾으면
  인터넷 없이 정확하고 재현된다. PubChem 은 (1) 원본에 지문이 없는 이름을 메우고
  (2) 사람이 읽을 대표 이름·CID 를 붙이는 보조로만 쓴다. 결과는 캐시에 남겨 다시 조회하지 않는다.

'왜' ChemSpider 는 안 쓰나?
  RSC API 키가 있어야 한다. 키가 생기면 pubchem_* 와 같은 모양의 함수를 덧대면 된다.

⚠️ 같은 지문 = 같은 구조 "후보"를 두 도구가 지목했다는 뜻이지 그 구조가 확정됐다는 뜻이 아니다.
"""
import gzip
import json
import re
import time
import csv
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path

CACHE = Path.home() / ".npdb_cache"          # build_candidate_db.py 와 같은 캐시(동기화 대상 아님)
IKEY_RE = re.compile(r"\b([A-Z]{14}-[A-Z]{10}-[A-Z])\b")
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


# ── MS-DIAL 사전(.msp) → 이름별 지문 ───────────────────────────────────────
def msp_index(msp_path) -> dict:
    """NAME → [InChIKey…]. 사전이 크므로(수십만 레코드) 한 번 읽어 캐시에 둔다.
    캐시 이름에 파일 크기·수정시각을 넣어 사전이 바뀌면 자동으로 다시 만든다."""
    msp = Path(msp_path)
    st = msp.stat()
    CACHE.mkdir(parents=True, exist_ok=True)
    cpath = CACHE / f"msp_ikey_{msp.stem}_{st.st_size}_{int(st.st_mtime)}.json.gz"
    if cpath.exists():
        with gzip.open(cpath, "rt", encoding="utf-8") as fh:
            return json.load(fh)
    idx, name, key = {}, None, None
    with msp.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("NAME:"):
                name, key = line[5:].strip(), None
            elif line.upper().startswith("INCHIKEY:"):
                key = line.split(":", 1)[1].strip() or None
            elif not line.strip():
                if name and key:
                    lst = idx.setdefault(name, [])
                    if key not in lst:
                        lst.append(key)
                name, key = None, None
    if name and key:
        idx.setdefault(name, [key])
    with gzip.open(cpath, "wt", encoding="utf-8") as fh:
        json.dump(idx, fh)
    return idx


# ── mzmine 후보 목록(LOTUS) → 이름별 지문 ─────────────────────────────────
def candidate_index(db_dir) -> dict:
    """name → [InChIKey…]. 같은 이름이 다른 구조를 가리킬 수 있어 목록으로 둔다."""
    idx = {}
    for f in sorted(Path(db_dir).glob("*_candidates_*.csv")):
        with f.open(encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                n, k = (r.get("name") or "").strip(), (r.get("inchikey") or "").strip()
                if n and k:
                    lst = idx.setdefault(n, [])
                    if k not in lst:
                        lst.append(k)
    return idx


# ── PubChem (보조) ────────────────────────────────────────────────────────
class PubChem:
    """이름→InChIKey, InChIKey→(CID, 대표 이름). 결과는 파일 캐시에 남긴다.
    '왜' 0.25 초 쉬나? PubChem 사용 규정이 초당 5회 이하다."""

    def __init__(self, cache_file=None, offline=False):
        self.path = Path(cache_file) if cache_file else CACHE / "pubchem_cache.json"
        self.offline = offline
        self.db = json.loads(self.path.read_text(encoding="utf-8")) if self.path.exists() else {}
        self.n_calls = 0

    def _get(self, url):
        if self.offline:
            return None
        time.sleep(0.25)
        self.n_calls += 1
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            return {} if e.code == 404 else None     # 404 = 없음(캐시함), 그 밖 = 일시 오류(캐시 안 함)
        except Exception:
            return None

    def _cached(self, key, url, parse):
        if key in self.db:
            return self.db[key]
        js = self._get(url)
        if js is None:
            return None
        val = parse(js) if js else None
        self.db[key] = val
        return val

    def ikey_from_name(self, name):
        q = urllib.parse.quote(name, safe="")
        return self._cached(f"name2ikey|{name}",
                            f"{PUBCHEM}/compound/name/{q}/property/InChIKey/JSON",
                            lambda js: js["PropertyTable"]["Properties"][0].get("InChIKey"))

    def title_from_ikey(self, ikey):
        return self._cached(f"ikey2title|{ikey}",
                            f"{PUBCHEM}/compound/inchikey/{ikey}/property/Title/JSON",
                            lambda js: {"cid": js["PropertyTable"]["Properties"][0].get("CID"),
                                        "title": js["PropertyTable"]["Properties"][0].get("Title")})

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.db, ensure_ascii=False, indent=0), encoding="utf-8")


# ── 이름 하나 → 지문 하나 ──────────────────────────────────────────────────
def ikey_for(name, index, pubchem=None):
    """(InChIKey 또는 None, 출처). 원본 사전에 하나로 정해지면 그것, 여럿이면 None('ambiguous'),
    없으면 이름 안에 박힌 지문, 그래도 없으면 PubChem."""
    if not name:
        return None, None
    keys = index.get(name, [])
    if len(keys) == 1:
        return keys[0], "source_db"
    if len(keys) > 1:
        blocks = {k[:14] for k in keys}
        return (keys[0], "source_db_stereo_variants") if len(blocks) == 1 else (None, "ambiguous_in_source")
    m = IKEY_RE.search(name)
    if m:
        return m.group(1), "in_name"
    if pubchem is not None:
        clean = name.split("!", 1)[-1].split(";", 1)[0].strip()   # "NCGC…!이름", "이름; CE30; …"
        k = pubchem.ikey_from_name(clean)
        if k:
            return k, "pubchem"
    return None, None


def relation(k1, k2, f1=None, f2=None, n1=None, n2=None):
    """두 도구 지목의 관계. 위에서부터 강한 일치."""
    norm = lambda s: "".join(ch for ch in (s or "").lower() if ch.isalnum())
    if k1 and k2:
        if k1 == k2:
            return "same_structure"
        if k1[:14] == k2[:14]:
            return "same_skeleton_stereo_differs"
    elif n1 and n2 and norm(n1) == norm(n2):
        return "same_name_no_key"
    if f1 and f2 and f1.replace(" ", "") == f2.replace(" ", ""):
        return "same_formula_different_structure"
    return "different_formula"
