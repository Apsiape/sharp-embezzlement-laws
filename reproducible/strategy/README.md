# Finite-chain Bell certificates

The scripts certify three finite XX-chain constructions. The shared resources and games are specified separately below. Run from this directory, with the dependencies in `../requirements.txt`, without Python's `-O` option. Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1` and unset `PYTHONOPTIMIZE` for the reference runs.

## Original game with one shared nine-dimensional register

```sh
python assisted_xx_certificate.py
python assisted_independent_check.py
```

The shared state is the centrally cut XX ground state tensored with one Φ₉. This register is an input resource, not a returned catalyst. The original game has 243 uniformly weighted question pairs. Providing the same register to a maximally entangled competitor leaves it maximally entangled.

| XX sites | Compressed covariance dimension | Certified interval for the prescribed strategy's loss |
| ---: | ---: | --- |
| 262,144 | 296 | [0.000746210786, 0.000747368515] |
| 524,288 | 320 | [0.000701624233, 0.000702783850] |
| 1,048,576 | 344 | [0.000661007046, 0.000662385080] |

The last two intervals lie below the universal flat-state loss floor `(3−2√2)/243 = 0.000706061215…`. The first is a control for this prescribed strategy; it does not exclude a better strategy at that length.

The primary verifier encloses exact Jacobi roots and quadrature weights, evaluates kernel entries by directed integer arithmetic, checks eigenvalue residuals using exact Rayleigh centers, and bounds the omitted singular mass. A geometric mode model is scored through exact integer rank counts. The three fixed rational amplitude-distance budgets are 0.00001363, 0.000009187, and 0.00001684. The program checks each budget and every error cap before using it to compute the loss interval. `assisted_xx_certificate.json` records the reference run; fields beginning with `nominal` are estimates, not certificates. Raw residual bounds can vary between numerical libraries; the fixed-budget loss intervals are portable when all checks pass.

The physical canonical strategy has loss `2 E(ψ_N)/243`. The square-root duplication-defect transfer proves the thresholds 0.000703 and 0.000663. The separately reported generic transfer applies the model POVMs to the physical state and has a different bound.

`assisted_independent_check.py` is a second implementation with a different kernel perturbation estimate, exact-rational mode ratios, 20 modes per parity, and a different Schmidt prefix. Both its canonical-strategy and fixed-POVM transfers independently certify losses below 0.000703 and 0.000663. Neither implementation expands the full physical Schmidt spectrum.

## Original game using the chain alone

```sh
python certify_xx84_modes.py
python model_score.py model84_exponents.json --cutoff 380000 --distance-bound 3/2500 --require-separation
```

The first script verifies covariance compression and the geometric model's vector distance from the actual open XX ground state at length `2^84`. Eigensolvers propose eigenpairs. The checks use exact rational residuals for the Jacobi matrices, binary64 arithmetic bounds for the compressed kernel, disjoint residual intervals for the selected eigenvalues, and a Frobenius-norm bound on omitted singular values. The paper derives the analytic error estimates.

The second script uses exact integer multiplicities, directed fixed-point arithmetic, and a classical strategy on the Schmidt tail. The model loss is below 0.000410720762. Transfer to the physical state gives a loss below 0.000460800. The `--require-separation` flag makes a failed strict separation exit unsuccessfully.

The covariance calculation has dimension 1418 and the model has 81 paired modes. The physical chain length is enormous; this result uses the chain as its sole shared resource.

## Four sites in the tagged game

```sh
python check_four_site_hybrid.py
```

SymPy checks the full spin Hamiltonian, normalized ground eigenvector, both local Schmidt-basis changes, and the Bell operator on the physical state. Exact radical bounds certify the score gaps. The four-site result uses a tagged game combining the original game with tilted CHSH. The strategy plays the original-game branch classically; the advantage comes from the tilted-CHSH branch. The chain ground state is the sole shared resource.
