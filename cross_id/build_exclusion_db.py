# -*- coding: utf-8 -*-
"""
S6c 준비 — '시료에서 나올 수 없는 구조'를 거르는 의약품 기록부 고정판 만들기 (2026-11-12, 이전 연구)
==========================================================================================
'왜' 필요한가?
  MSNovelist(SIRIUS 6 의 새 구조 생성)는 목록에 없는 구조를 만들어 낸다. 만든 구조가 이미 알려진
  의약품과 같거나 의약품 쪽에 더 닮았으면, 식물 추출물에서 나왔다고 보기 어렵다. 그 판단을
  사람의 기억이 아니라 공개 기록부와의 대조로 남기기 위해 기록부를 날짜 박힌 고정판으로 받아 둔다.

기록부
  ① DrugCentral (기본) — 미국·유럽·일본 승인 의약품과 동물용 의약품의 구조 파일(SDF, MOL V2000).
     내려받기 쪽: https://drugcentral.org/download  (판 = 주소의 연도 폴더, 2026-11-12 기준 '2027')
     약 4,300 구조. 천연물 유래 약물(모르핀·파클리탁셀 등)도 섞여 있으므로 판정은 하류에서
     LOTUS 기록과 나란히 본다(LOTUS 에도 있으면 시료 유래 가능).
  ② ChEMBL 임상 단계 이상(max_phase ≥ 1, --chembl) — 개발 중 약물까지 넓힌다. 천연물 표지(natural_product)가 있다.
     2026-11-12 밤 ChEMBL 웹 API 가 오류(500)를 내 받지 못했다. 서버가 살아나면 --chembl 로 추가한다.

'왜' 고정판인가? 기록부는 판이 바뀐다. 같은 후보가 다음 달 다른 판정을 받지 않도록 받은 날짜·판을 남긴다.

산출: db_frozen/drugs_<출처>_<판>_<받은날>.csv.gz  +  같은 이름 .json
      칸: source, drug_id, name, max_phase, natural_product, smiles, inchikey
쓰는 법: py scripts/build_exclusion_db.py [--chembl]      (인터넷 필요, RDKit 필요)
"""
import sys, io, json, gzip, csv, time, datetime, argparse, urllib.request, urllib.parse
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "db_frozen"            # 프로젝트 안(수 MB) — 다른 PC 로 함께 넘어가야 판정이 같아진다
DC_RELEASE = "2027"
DC_URL = f"https://unmtid-dbs.net/download/DrugCentral/{DC_RELEASE}/structures.molV2.sdf.gz"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
UA = {"User-Agent": "natural-product-comparison-pipeline/1.0"}   # 개인정보 없는 일반 이름(전역 규칙)
COLS = ["source", "drug_id", "name", "max_phase", "natural_product", "smiles", "inchikey"]


def fetch(url, timeout=180):
    """주소 하나를 받는다. 일시 오류는 세 번까지 다시 시도한다."""
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return r.read()
        except Exception as e:
            if i == 2:
                raise
            print(f"  다시 시도 {i+1}: {e}")
            time.sleep(5 * (i + 1))


def save(rows, stem, meta):
    OUT_DIR.mkdir(exist_ok=True)
    with gzip.open(OUT_DIR / f"{stem}.csv.gz", "wt", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLS)
        w.writerows(rows)
    meta.update({"n_structures": len(rows), "columns": COLS})
    (OUT_DIR / f"{stem}.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"저장: {OUT_DIR / (stem + '.csv.gz')}  ({len(rows)} 구조)")


def drugcentral(today):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    raw = gzip.decompress(fetch(DC_URL))
    sup = Chem.ForwardSDMolSupplier(io.BytesIO(raw))
    rows, n_bad = [], 0
    for m in sup:
        if m is None:                      # 읽히지 않는 구조(금속 착물 등)는 세어 두고 건너뛴다
            n_bad += 1
            continue
        rows.append(["DrugCentral", m.GetProp("ID") if m.HasProp("ID") else "",
                     m.GetProp("PREFERRED_NAME") if m.HasProp("PREFERRED_NAME") else "",
                     "", "", Chem.MolToSmiles(m), Chem.MolToInchiKey(m)])
    save(rows, f"drugs_DrugCentral_{DC_RELEASE}_{today}",
         {"source": "DrugCentral structures.molV2.sdf.gz", "url": DC_URL, "release": DC_RELEASE,
          "downloaded": today, "n_unreadable_skipped": n_bad,
          "license_note": "drugcentral.org/download 의 License 절 확인 필요(2026-11-12 미확인)"})


def chembl(today):
    status = json.loads(fetch(f"{CHEMBL_API}/status.json"))
    release = status.get("chembl_db_version", "unknown")
    url = f"{CHEMBL_API}/molecule.json?" + urllib.parse.urlencode({"max_phase__gte": 1, "limit": 1000, "format": "json"})
    rows = []
    while url:
        d = json.loads(fetch(url))
        for m in d["molecules"]:
            st = m.get("molecule_structures") or {}
            if st.get("canonical_smiles") and st.get("standard_inchi_key"):
                rows.append(["ChEMBL", m["molecule_chembl_id"], m.get("pref_name") or "", m.get("max_phase"),
                             m.get("natural_product"), st["canonical_smiles"], st["standard_inchi_key"]])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
        print(f"  누적 {len(rows)}")
    save(rows, f"drugs_ChEMBLclinical_{release}_{today}",
         {"source": "ChEMBL web services molecule?max_phase__gte=1", "release": release, "downloaded": today,
          "license": "CC BY-SA 3.0"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chembl", action="store_true", help="ChEMBL 임상 단계 이상 약물도 받는다")
    args = ap.parse_args()
    today = datetime.date.today().isoformat()
    drugcentral(today)
    if args.chembl:
        chembl(today)


if __name__ == "__main__":
    main()
