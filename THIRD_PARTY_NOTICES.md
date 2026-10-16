# Source and template attribution

The files in `manuscript/` named `optica-article.cls`, `styles/opticajournal.sty`, `opticajnl.bst`, `jabbrv.sty`, and the two language files were copied from the official Optica universal template archive:

https://opg.optica.org/resources/author/templates/latex-universal-template.zip

The `manuscript/osa-supplemental-document.cls` was copied from:

https://opg.optica.org/resources/author/templates/latex-supplemental-document.zip

Downloaded 16 October 2026. Their original notices are retained. The class files identify the LaTeX Project Public License version 1.3 or later. Other bundled files retain their own headers. Inclusion does not imply endorsement of the research by Optica or Overleaf.

The exact fixture incorporates tabulated CIE data from the International Commission on Illumination:

- CIE (2019), *Colour-matching functions of CIE 1931 standard colorimetric observer*, DOI 10.25039/CIE.DS.xvudnb9b.
- CIE (2019), *CIE standard illuminant D65*, DOI 10.25039/CIE.DS.hjfjmt59.
- CIE (2018), *CIE standard illuminant A - 1 nm*, DOI 10.25039/CIE.DS.8jsxjrsn.
- CIE (2018), *Relative spectral power distributions of illuminants representing typical fluorescent lamps, 1 nm wavelength steps*, DOI 10.25039/CIE.DS.54hy6srn.

The retained source metadata identify [Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/). The fixture selects 380 to 780 nm and the FL11 column, converts finite decimals to rational values, uses piecewise-linear interpolation, and normalizes each illuminant to continuous white Y=100. These transformations define the research model and are not CIE endorsements or new measurements. Preserve attribution and applicable share-alike obligations for adapted data in any later distribution. This notice does not assign a license to the entire original manuscript or project code.

The portable compiler used for local build verification was Tectonic 0.17.0 from its official GitHub release. Its downloaded Windows MSVC archive matched the publisher's SHA-256 digest `f61ce51f0b0ade1015b7de7ef368541c5424e9756ecbd0d7af97d6d48030845f`. Compiler binaries and package caches are excluded from the manuscript and Git backup. Manuscript sources were compiled locally; no Overleaf upload was performed.
