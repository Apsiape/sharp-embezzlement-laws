# Sharp embezzlement laws for quantum chains and finite Bell tests

Seth Douglas and Nidhal Mghirbi

Version 1.0.0, November 2026. The complete paper is `sharp-embezzlement-laws.pdf`. LaTeX source and all four figures are in `manuscript/`. Programs, datasets, and arithmetic certificates are in `reproducible/`. `CHANGES.md` records this version and revision 4, `REVIEW_RESPONSE.md` explains the treatment of the 11 November 2026 review feedback, and `VERIFICATION.md` records the checks and their scope.

## Main results

- Sharp optimal local-unitary embezzlement laws for finite open XX and Ising chains, with an inverse-square-logarithm remainder, arbitrary open cuts, and a growing-target Gaussian crossover. The crossover curve is credited to earlier LOCC work; the chain theorem controls the exact discrete spectra and local-unitary optimization.
- An explicit loss floor `(3−2√2)/243 ≈ 0.000706061215` for maximally entangled states of every finite dimension in the uniform 243-pair game. The qualitative floor phenomenon is prior work.
- With one shared nine-dimensional maximally entangled register, the open XX ground state at **524,288 sites** admits loss below **0.000703**. At **1,048,576 sites**, its certified loss is below **0.000663**.
- With the chain as the sole shared resource, separation is certified at length `2^84`, with loss below **0.000461**. The sharper computed endpoint remains in the certificate. The four-spin separation realizes a known tilted-CHSH witness in a separately defined tagged game, with classical play in its original-game branch.
- A dimension-independent lower comparison with the explicit constant **1,353,145.5**, covering arbitrary local POVMs.
- A sharp spectral law for the **imported XXZ corner-transfer-matrix spectrum** as the anisotropy approaches one from above. The spectral formula's physical identification is an input. Exact circuit conjugates of XX provide separate interacting bulk examples.

## Build the paper

With a TeX Live installation containing the source's standard packages:

```sh
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error type-of-matter.tex
```

Alternatively, run `pdflatex type-of-matter.tex` three times. The bibliography is included in the main source.

## Verify the package and calculations

From the archive root:

```sh
python verify_manifest.py
python -m pip install -r reproducible/requirements.txt
sh verify.sh
```

`verify.sh` checks every manifest entry, the nine regressions, the algebraic certificates, both assisted-chain implementations, the chain-only certificate, and the exact rounding and threshold inequalities. It exits unsuccessfully if a required comparison fails. The runner disables Python optimization and fixes the BLAS thread count. SymPy is a declared dependency for the four-spin certificate.

The source manifest covers every distributed payload file, including the PDF and reference verification output. It excludes only itself and the repository metadata files `CITATION.cff` and `.zenodo.json`. A build or regenerated dataset may add or change files; the manifest identifies the distributed bytes.

Figure and numerical-data regeneration commands are in `reproducible/README.md`; certificate parameters and commands are in `reproducible/strategy/README.md`. The records in `verification/` describe the AI-assisted proof audits and their scope. This release contains no Lean formalization.

## Licence

Text and documentation: CC BY 4.0. Code: MIT or Apache 2.0, at your option. See `RIGHTS.md`.
