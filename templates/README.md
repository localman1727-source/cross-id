# mzmine batch templates

`example_waters_mse_POS.mzbatch` / `example_waters_mse_NEG.mzbatch` are the peak-picking templates CROSS-ID starts from
(mzmine 4; Waters SYNAPT XS, MSe acquisition). `make_mzmine_batches.py` reads one of them, swaps in your raw files,
polarity and output paths, and adapts the steps for phytochemical profiling.

They are processing settings only — no raw data. Tune mass-detection noise levels, RT range and alignment
tolerances to your own instrument and method, then save your own copy here (e.g. `my_method_POS.mzbatch`)
and point `processing.mzmine_template` in the study config to it.
