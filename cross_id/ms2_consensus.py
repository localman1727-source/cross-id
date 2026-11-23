# -*- coding: utf-8 -*-
"""
봉우리별 조각무늬(MS2) 합본 만들기 — 여러 런에서 반복된 조각만 남긴다
=====================================================================
'왜' 합본인가?
  mzmine 의 sirius.mgf 에는 봉우리마다 런별 MS2 가 최대 12개 있다. 한 런만 쓰면 그 순간 섞인 잡음
  봉우리가 들어가고 약한 진짜 조각이 빠질 수 있다. 진짜 조각은 여러 런에서 같은 자리에 반복되므로,
  반복된 것만 남기면 더 깨끗한 무늬가 된다.

'왜' 합치기 전에 닮음을 보나?
  런마다 무늬가 크게 다르면(물질 두 개가 겹쳤거나 특정 시료에만 다른 물질이 섞임) 그대로 합칠 때
  두 물질의 조각이 섞인 가짜 무늬가 된다. 그래서
    1) 각 런 무늬를 가장 큰 조각 = 1 로 맞추고
    2) 런끼리 코사인 닮음(같은 질량 ±MZ_TOL 조각끼리 짝지어 계산)을 구해
    3) 서로 닮은(코사인 ≥ OUTLIER_COS) 런이 가장 많이 모인 무리 하나만 남기고
       ('왜' 전체 중앙값이 아닌가? 10-08 실측: 어떤 산지 시료에서 물질이 약하면 그 런의 무늬는 조각 2~3개짜리
        잡음이고, 근처에 다른 물질이 겹쳐 들어온 런도 있다. 이런 런이 절반 가까이면 '전체와의 중앙값'
        기준으로는 멀쩡한 런까지 모두 탈락한다 — 1,994개 중 585개가 그렇게 합본을 못 만들었다.)
    4) 남은 런에서 절반 이상(MIN_FRACTION)에 나타난 조각만 남긴다(질량 = 세기 가중 평균, 세기 = 평균).
  닮은 짝이 하나도 없으면 합본 대신 신호가 가장 센 런 하나를 쓰고 표시한다(mode 칸).
  기준값(0.6, 절반)은 이 자료에서 정한 것이 아니라 흔히 쓰는 출발값이다 — 결과는 한 런 방식과
  나란히 비교해 판단한다(run_ms2query.py --spectrum last / consensus).

순수 파이썬(외부 패키지 없음) — MS-FINDER 입력 스크립트(py 3.13)와 MS2Query venv(3.11) 양쪽에서 쓴다.
"""
import math
from pathlib import Path

MZ_TOL = 0.01          # Da. 같은 조각으로 볼 질량 차(Waters QTof MS2)
OUTLIER_COS = 0.6      # 이 이상이면 '서로 닮은 런'으로 본다
MIN_FRACTION = 0.5     # 남은 런 중 이 비율 이상에 나타난 조각만 남긴다


def read_ms2_blocks(mgf):
    """FEATURE_ID → [{pmz, rt_sec, file, adduct, peaks[(mz, int)]}...] (런별 MS2 전부, 파일 순서)."""
    out = {}
    text = Path(mgf).read_text(encoding="utf-8", errors="replace")
    for blk in text.split("END IONS"):
        lines = [l.strip() for l in blk.strip().splitlines() if l.strip()]
        if not lines:
            continue
        head = dict(l.split("=", 1) for l in lines if "=" in l and not l[:1].isdigit())
        if head.get("MSLEVEL") != "2":
            continue
        pmz = float(head["PEPMASS"].split()[0])
        peaks = []
        for l in lines:
            if l[:1].isdigit():
                a = l.split()
                m, i = float(a[0]), float(a[1])
                if m <= pmz + 0.5 and i > 0:          # 전구 이온 위쪽 잡음 제거
                    peaks.append((m, i))
        if peaks:
            out.setdefault(head["FEATURE_ID"], []).append({
                "pmz": pmz, "rt_sec": float(head["RTINSECONDS"]), "file": head.get("FILENAME", ""),
                "adduct": head.get("ADDUCT", ""), "peaks": peaks})
    return out


def _norm(peaks):
    mx = max(i for _, i in peaks)
    return sorted((m, i / mx) for m, i in peaks)


def cosine(a, b, tol=MZ_TOL):
    """정규화한 두 무늬의 코사인 닮음. 큰 조각부터 질량 ±tol 안에서 1:1 로 짝짓는다."""
    pairs = []
    j0 = 0
    for ia, (ma, xa) in enumerate(a):
        while j0 < len(b) and b[j0][0] < ma - tol:
            j0 += 1
        j = j0
        while j < len(b) and b[j][0] <= ma + tol:
            pairs.append((xa * b[j][1], ia, j))
            j += 1
    used_a, used_b, s = set(), set(), 0.0
    for p, ia, jb in sorted(pairs, reverse=True):
        if ia not in used_a and jb not in used_b:
            used_a.add(ia); used_b.add(jb); s += p
    na = math.sqrt(sum(x * x for _, x in a)); nb = math.sqrt(sum(x * x for _, x in b))
    return s / (na * nb) if na and nb else 0.0


def _median(v):
    v = sorted(v)
    n = len(v)
    return (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2) if n else float("nan")


def consensus(spectra):
    """런별 MS2 목록 → (합본 peaks, 기록 dict)."""
    normed = [_norm(s["peaks"]) for s in spectra]
    n = len(normed)
    if n == 1:
        return normed[0], {"n_spectra": 1, "n_used": 1, "median_cosine": None, "mode": "single_only",
                           "used_files": spectra[0]["file"]}
    cos = [[1.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            cos[i][j] = cos[j][i] = cosine(normed[i], normed[j])
    tic = [sum(i for _, i in s["peaks"]) for s in spectra]          # 원래 세기 합(정규화 전) = 신호 세기
    nbrs = [{j for j in range(n) if j != i and cos[i][j] >= OUTLIER_COS} for i in range(n)]
    center = max(range(n), key=lambda i: (len(nbrs[i]), tic[i]))   # 닮은 런이 가장 많은 런(동률이면 신호 센 쪽)
    keep = sorted({center} | nbrs[center])
    info = {"n_spectra": n, "n_used": len(keep),
            "median_cosine": round(_median([cos[i][j] for i in range(n) for j in range(i + 1, n)]), 4),
            "used_files": ";".join(spectra[k]["file"] for k in keep)}
    if len(keep) < 2:
        best = max(range(n), key=lambda i: tic[i])
        info["mode"] = "inconsistent_fallback_single"
        info["used_files"] = spectra[best]["file"]
        return normed[best], info
    # 남은 런의 조각을 질량순으로 늘어놓고 ±tol 안에서 이어지는 것끼리 한 묶음
    allp = sorted((m, x, k) for k in keep for m, x in normed[k])
    groups, cur = [], [allp[0]]
    for p in allp[1:]:
        if p[0] - cur[-1][0] <= MZ_TOL:
            cur.append(p)
        else:
            groups.append(cur); cur = [p]
    groups.append(cur)
    need = max(2, math.ceil(MIN_FRACTION * len(keep)))
    out = []
    for g in groups:
        runs = {k for _, _, k in g}
        if len(runs) >= need:
            w = sum(x for _, x, _ in g)
            out.append((sum(m * x for m, x, _ in g) / w, w / len(keep)))
    if not out:                                         # 반복 조각이 하나도 없으면 가장 조각 많은 런
        best = max(keep, key=lambda i: tic[i])
        info["mode"] = "no_repeated_peaks_fallback_single"
        info["used_files"] = spectra[best]["file"]
        return normed[best], info
    info["mode"] = "consensus"
    return _norm(out), info
