"""설치된 질량분석 도구가 실제 계산과 파일 왕복을 수행하는지 확인한다."""
import sys
import io

# Windows에서도 한국어 기록이 깨지지 않도록 출력 인코딩을 지정한다.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import json
import tempfile
from pathlib import Path
from importlib.metadata import version

import numpy as np
import pyopenms as oms
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from matchms import Spectrum
from matchms.similarity import CosineGreedy


def main():
    """가상 신호로 계산을 확인하고 실제 연구자료에는 접근하지 않는다."""
    result = {"python": sys.executable,
              "versions": {name: version(name) for name in
                           ("matchms", "pyopenms", "rdkit", "numpy", "scipy")}}

    # RDKit이 구조를 읽고 올바른 분자식을 만드는지 확인한다.
    molecule = Chem.MolFromSmiles("CCO")
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    assert formula == "C2H6O", formula
    result["rdkit_formula"] = formula

    # 동일한 조각 신호의 대조 결과가 완전 일치하는지 확인한다.
    masses = np.array([100., 150., 200.])
    intensities = np.array([0.2, 0.5, 1.0])
    spectrum = Spectrum(mz=masses, intensities=intensities,
                        metadata={"precursor_mz": 250.0})
    match = CosineGreedy(tolerance=0.01).pair(spectrum, spectrum)
    score = float(match["score"])
    count = int(match["matches"])
    assert np.isclose(score, 1.0) and count == 3, str(match)
    result["matchms_self_comparison"] = {"score": score, "matches": count}

    # OpenMS의 네이티브 파일 입출력까지 시험해 import만 되는 상태와 구분한다.
    sample = oms.MSSpectrum()
    sample.setMSLevel(1)
    sample.set_peaks((masses, intensities))
    experiment = oms.MSExperiment()
    experiment.addSpectrum(sample)
    # OpenMS Windows 파일 입출력에는 한글이 없는 임시 경로를 사용한다.
    temp_root = "C:/massspec_env/tmp" if sys.platform == "win32" else None
    with tempfile.TemporaryDirectory(prefix="massspec_verify_", dir=temp_root) as folder:
        filename = str(Path(folder) / "synthetic.mzML")
        oms.MzMLFile().store(filename, experiment)
        loaded = oms.MSExperiment()
        oms.MzMLFile().load(filename, loaded)
        actual_mz, actual_i = loaded[0].get_peaks()
        assert len(loaded) == 1 and np.allclose(actual_mz, masses)
        assert np.allclose(actual_i, intensities)
    result["pyopenms_mzml_roundtrip"] = "passed"
    result["status"] = "passed"

    # 다른 세션에서도 검증 근거를 확인할 수 있도록 프로젝트 안에 남긴다.
    target = Path(__file__).resolve().with_name("verify_massspec_environment_results.json")
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
