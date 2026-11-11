# Numerical reproduction and exact certificates

This directory accompanies *Sharp embezzlement laws for quantum chains and finite Bell tests*. It contains the programs and data for the chain, conformal, and repeated-pair figures, the two finite-state witnesses for the 243-pair game, and the exact arithmetic certificates. The programs run with Python 3.10 or later and the packages in `requirements.txt`. The release checks used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, and SymPy 1.14.0.

## Run the checks and regenerate the figures

Run these commands from this directory.

```sh
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONOPTIMIZE
python test_reproducibility.py
python certificate.py
python constant_verify.py
python compare_data.py
python plot.py
```

`plot.py` writes `fig_main.pdf`, `fig_law.pdf`, `fig_diagnostics.pdf`, `figure_data.csv`, and two TeX table fragments to `generated/figures`. `python plot.py --archive-only --output archived-figures` regenerates the first two plots and the table fragments from the archived data. The bundled figure data use regenerated critical and conformal results and the archived gapped sweeps. The diagnostics plot reports the loss of the explicit measurements. The dashed curves in its right panel use the shared state Φ₉.

`certificate.py` verifies the 19 Schur-pivot intervals for the strategy upper bound \(7\mathcal E/3\). `constant_verify.py` checks the rational Gram certificate and scalar comparisons for the flagged-sector coefficient \(210\) and the lower-comparison constant \(1{,}353{,}145.5\).

## Finite-state witnesses

The uniform loss floor for maximally entangled states is \((3-2\sqrt2)/243 \simeq 7.06061215\times10^{-4}\). Both witnesses use a shared Φ₉ in the specified state, with the ranks below including that factor.

| Witness | Schmidt rank | Qubits per party | Explicit strategy loss | Data |
| --- | ---: | ---: | ---: | --- |
| Six integer shell coefficients | 288 | 9 | \(6.88936502\times10^{-4}\) | `generated/witness-compact.json` |
| Dyadic sine shells, \(K=9\) | 4,608 | 13 | \(4.02257840\times10^{-4}\) | `generated/witness.json` |

`python witness.py --compact` uses exact \(\mathbb Q(\sqrt2)\) arithmetic to verify the sorted spectrum, shell norm \(555\), duplication error \((598-390\sqrt2)/555\), and strict loss gap. `python witness.py` checks the sine shell sum against compressed rank intervals and reports a proved loss upper bound of \(5.64106333\times10^{-4}\). The upper bound is \(6.63110024\times10^{-4}\) after allowing a phase-aligned vector preparation error of \(0.002\).

## Regenerate numerical data

The following commands regenerate the bundled data. The full critical sweep evaluates each Ising size once and doubles its one-body spectrum to obtain the corresponding XX point. On the recorded machine the final Ising/XX pair took 105 seconds and the preceding pair took 134 seconds. FFT and Lanczos timing varies with system size and the installed numerical libraries.

```sh
python run.py chains critical --min-power 3 --max-power 17 --output generated/critical.json
python run.py theory --output generated/theory.json
python run.py pairs --output generated/pairs.json
python run.py gapped-floor 0.5 0.1 0.01 --output generated/gapped-floor-samples.json
python witness.py > generated/witness.json
python witness.py --compact > generated/witness-compact.json
python compare_data.py > generated/archive_comparison.json
python plot.py
```

For a shorter computation, `python run.py chains critical --max-power 10 --output small-critical.json` generates both critical models through Ising size 1024 and XX size 2048. `python run.py chains xx --max-power 10 --output direct-xx.json` computes the XX points directly and provides an independent check of the doubled modes. The `--modes` and `--threshold` options permit sensitivity checks.

These commands regenerate the original gapped sweep sizes and floor estimates.

```sh
python run.py chains gapped --delta 0.1 --max-power 16 --output generated/gapped-0.1.json
python run.py chains gapped --delta 0.01 --max-power 16 --output generated/gapped-0.01.json
python run.py gapped-floor 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002 --output generated/gapped-floor.json
python plot.py
```

The floor runner compares successive chain lengths and records whether the change in \(\mathcal E\) on doubling falls below \(10^{-6}\). It also records when the maximum size was reached. The bundled regeneration ran the three floor samples above. The original large-size gapped points remain available in the archived files.

## Calculations and error accounting

- `chains.py` computes dense or matrix-free one-body spectra. For large chains it uses a discrete sine transform and Lanczos eigenvalues of a positive matrix with eigenvalues that identify the entangled modes. `spectrum.py` enumerates products above relative threshold \(10^{-11}\). `E` is the homogeneous squared defect of the retained probabilities and matches the original table convention. `E_normalized` divides it by their mass. The original program named `H_functional` returns twice the squared defect.
- `E_lower` and `E_upper` use the normalized Schmidt-tail inequality with the mass missing from the computed mode product. `E_lower_with_modes` and `E_upper_with_modes` also allow for omitted one-body modes. This last estimate treats the smallest computed occupation as an upper bound on each omitted occupation and assumes the retained eigenvalues are the leading exact ones. It measures spectral truncation under that assumption. It does not certify Lanczos eigenvalue or floating-point errors. The tests compare a sample of dense, matrix-free, doubled, and exact spin spectra.
- `loss_explicit` is the exact-rule loss formula evaluated on the normalized retained spectrum. `loss_B9` is the Hellinger block-averaging upper bound evaluated on the same spectrum. `loss_upper_prefix` uses the better of the explicit strategy and the classical strategy on the retained Schmidt support, then the classical strategy on the omitted support. Its missing-mode qualification is the same as for the one-body tail estimate. Neither the two probabilities nor the prefix bound is a numerical optimization over all strategies.
- `conformal.py` counts exact discrete tower degeneracies using integer arithmetic and computes the defect by intersecting rank intervals. There is no expansion of the individual Schmidt list. The tower data extend to \(\varepsilon=0.02\) for each tower. An exponential-moment bound controls the discarded probability above the retained energy cutoff. The partition function and entropy have independent geometric-series tail bounds. The smooth curve uses Bessel inversion and adaptive quadrature in logarithmic rank, with an analytic right-tail mass bound and a separate quadrature estimate. These tail bounds do not include floating-point rounding.
- The pair calculation also intersects exact integer rank intervals, with exact binomial degeneracies.
- `game.py` constructs the Pauli, sign, and label projectors and directly evaluates all 243 question pairs on small spectra. The tests compare that calculation to the quadratic score formula, verify the classical score, and recompute the nineteen Schur pivots in exact \(\mathbb{Q}(\sqrt2)\) arithmetic.

`original/supplement` preserves the files from the original supplement ZIP with SHA-256 `8eebb2853fb0eb26988e3e985ef103d67fb699a16c5e3870edd9effe5c8d58a8`. Its exploratory scripts retain their original machine-specific paths and references to intermediate files outside the archive. Use the portable scripts above to read the archived JSON data or recompute the intermediates from model parameters.


## Finite critical-chain certificates

The `strategy` directory contains certificates for the assisted and chain-only results in the original game, and for the four-spin result in the tagged game. Run without Python's `-O` option:

```sh
python strategy/check_four_site_hybrid.py
OPENBLAS_NUM_THREADS=1 python strategy/assisted_xx_certificate.py
OPENBLAS_NUM_THREADS=1 python strategy/assisted_independent_check.py
OPENBLAS_NUM_THREADS=1 python strategy/certify_xx84_modes.py
python strategy/model_score.py strategy/model84_exponents.json --cutoff 380000 --distance-bound 3/2500 --require-separation
```

With one shared nine-dimensional maximally entangled register, the original game is certified on the exact XX ground state at 524,288 sites with loss below 0.000703, and at 1,048,576 sites with loss below 0.000663. The primary verifier checks covariance compression, quadrature roots and weights, integer eigenpair residuals, spectral truncation and Schmidt-tail errors. It uses fixed rational distance budgets to obtain portable outward score intervals. Its JSON separates certified bounds from nominal floating-point estimates. The 262,144-site control bounds the prescribed strategy's loss, not the optimal game loss. The independent verifier uses a different exact-rational mode model and prefix, and checks both the canonical physical-state strategy and transfer of the model POVMs.

The four-spin certificate uses exact symbolic arithmetic on the physical spin Hamiltonian and Bell observables. It concerns the tagged tilted-CHSH game, with classical play in its original-game branch. The original uniform game is certified on the bare XX ground state at length \(2^{84}\): loss below \(0.000460800\), compared with the universal maximally-entangled-state floor above \(0.000706061\). The physical length is enormous; the covariance verifier has dimension 1418. Its eigenvalue residual bounds, binary64 rounding bounds, and analytic compression error are derived in the paper. The score verifier uses integer multiplicities and directed arithmetic. The proof and the complete list of model exponents are included in the paper.

## Interacting XXZ data

```sh
python xxz/xxz_numerics.py > xxz/xxz_numerics.csv
python xxz/plot_xxz.py --output generated/figures/fig_xxz.pdf
```

The 21 data rows use exact integer level multiplicities for seven positive anisotropy spacings and three target dimensions. The omitted-tail bound is below \(10^{-15}\) at every point. Gibbs weights and plotted values use floating-point arithmetic; this plot is a numerical check of the analytic theorem. The limit is \(\Delta\downarrow1\) through pure gapped Néel phases.
