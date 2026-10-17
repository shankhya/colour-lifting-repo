# Response Tolerance and Spectral Regularity in Reference-Matched Reflectance Design
Shankhya Debnath
Department of Printing Technology, Regional Institute of Printing Technology, Kolkata, India.

Manuscript, supplement and reproducibility files.

## Files

- [`manuscript/`](manuscript/): `main.tex` and `main.pdf`, `supplement.tex` and `supplement.pdf`, bibliography, four vector figures, and the required Optica template files.
- [`reproducibility/`](reproducibility/): exact input data and certificates, Python verification scripts, MATLAB plotting code, three numerical tables, and the validation record.
- [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md): CIE data attribution and template license notices.
- [`SHA256SUMS`](SHA256SUMS): file integrity checksums.

## Compile the manuscript

Upload the **contents of `manuscript/`**, preserving its `figures/` and `styles/` subfolders, to an Overleaf project. Select **XeLaTeX** and **`main.tex`** as the main document. To build the supplement, select **`supplement.tex`** instead. Both entry points use the same local template files. A separate `main-review.tex` is unnecessary.

For a local TeX installation, run from `manuscript/`:

```text
latexmk -xelatex main.tex
latexmk -xelatex supplement.tex
```

The article has 10 pages and the supplement has 8 pages. The PDFs are included for reading without a TeX installation. Bibliography output and ZIP archives are generated as needed rather than stored as duplicate source packages.

## Verify the results

From the repository root, using Python 3.10 or later:

```text
python -B reproducibility/tools/verify_package.py
python -B reproducibility/tools/verify_paper_results.py
```

Both commands are read-only and use only Python's standard library. They require no network, optimizer or original research repository. Do not use `python -O`; the scientific checker requires assertions. The exact check takes about 20 seconds on the preparation machine.

Successful verification reports an absolute-bonus ratio lower bound of about `1.5002782835820574` and a certified tolerance interval `[0.056736, 0.137809]`. It reconstructs the integrated inputs and checks rank, exact feasible spectra, all-mixture constraints, full-class supports and the quantitative comparison/response certificates. The smooth-spectrum and input-continuity proofs are in the supplement.

The frozen fixture is `reproducibility/results/research/R135-paper-fixture.json`; the saved numerical export is `R135-paper-verified.json` in the same directory. The `R135` identifiers preserve source traceability.

## Reproduce figures and tables

In MATLAB, change to `reproducibility/` and run:

```matlab
run_paper_results('verify')
run_paper_results('render')
```

The code uses base MATLAB and a JVM for SHA-256; no optimization toolbox is required. The supplied figures were rendered with MATLAB R2025b. `verify` checks the fixture hash and saved numerical export without writing files. `render` writes four PDF/PNG figures and three CSV tables under `reproducibility/results/figures/paper/`. The distributed copies are in `manuscript/figures/` and `reproducibility/tables/`.

To regenerate the numerical export from the exact checks, run from `reproducibility/`:

```text
python -B tools/verify_paper_results.py --emit results/research/R135-paper-verified.json
```

This optional command overwrites the saved export and records a new elapsed time, changing its checksum. Run the read-only verification commands first. Figure binary hashes can vary with rendering software; the certified numerical quantities are the reproducibility target.
