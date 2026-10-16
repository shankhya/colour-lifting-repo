# Reference-matched reflectance design: reproducibility package

**Response Tolerance and Spectral Regularity in Reference-Matched Reflectance Design**

This private repository contains the revised paper, detailed proofs, exact computational certificates, and figure code. It is a standalone extraction of the research project, with a fresh distribution history. The manuscript remains an author-review draft; author names, affiliations and declarations are blank. No journal submission or reviewer access has been arranged.

## Contents

- `main.pdf`, `main.tex`, `main.bbl`, `references.bib`: revised article and bibliography.
- `main-review.tex`, `main-review.pdf`: matching Overleaf entry point and PDF.
- `overleaf/colour-lifting-josaa-overleaf.zip`: minimal 15-file journal source package. Select XeLaTeX and compile `main-review.tex` and `supplement.tex` separately.
- `supplement.pdf`, `supplement.tex`: detailed analytical proofs and certificate definitions.
- `figures/`, `tables/`: four vector figures and three numerical tables used by the paper.
- `reproducibility/`: exact CIE-derived input fixture, feasible spectra, support certificates, saved display values, Python verifier and MATLAB plotting function.
- `THIRD_PARTY_NOTICES.md`: source attribution and license terms.
- `validation.json`, `SHA256SUMS`, `verify_package.py`: distribution validation and integrity checks.

Internal research logs, schedules, credentials, machine configuration, caches, unrelated experiments and the original Git history are outside this distribution. Scientific inputs and unsuccessful-solver provenance inside the exact fixture are preserved.

## Reproduce the certified results

From the repository root, with Python 3.10 or later (tested with 3.12.14):

```text
python -B verify_package.py
python -B reproducibility/tools/verify_paper_results.py
```

Both commands are read-only. They need only Python's standard library, no network, optimizer, original repository or personal configuration. Do not use `python -O`; the scientific checker refuses disabled assertions. The checker takes about 20 seconds on the preparation machine; runtime varies.

Successful verification reports an absolute-bonus ratio lower bound of about `1.5002782835820574` and a certified tolerance interval `[0.056736, 0.137809]`. It reconstructs the integrated inputs and checks rank, exact feasible spectra, all-mixture constraints, full-class supports, and the quantitative comparison/response certificates. The smooth-spectrum and input-continuity results are analytical proofs in the supplement. No new optimizer search is needed to verify the certified bounds.

The frozen input fixture is `reproducibility/results/research/R135-paper-fixture.json`. Its source attribution and interpolation/normalization conventions are in the supplement and notices. The `R135` identifier is retained for exact source traceability.

## Reproduce plots and tables

In MATLAB, change to `reproducibility/`, then run:

```matlab
run_paper_results('verify')
run_paper_results('render')
```

The saved figures were rendered with MATLAB R2025b. The code uses base MATLAB, a JVM for SHA-256, and no optimization toolbox. `verify` checks the saved presentation export and fixture hash without writing files. `render` writes four PDF/PNG plots and three tables under `reproducibility/results/figures/paper/`, leaving the distributed figures unchanged. These are plots of mathematical certificates, not measurements. Figure binary hashes can depend on rendering software; the certified numerical quantities are the reproducibility target.

To explicitly regenerate the display JSON from exact checks, run from `reproducibility/`:

```text
python -B tools/verify_paper_results.py --emit results/research/R135-paper-verified.json
```

This optional command overwrites the display export and records a new elapsed time, so its checksum will differ. Start with the read-only commands above to validate the distributed files.

The minimal Overleaf ZIP contains only the two TeX sources, bibliography, four vector figures and required template dependencies. The reproducibility files are available separately in this repository. The title update changes no scientific result.
