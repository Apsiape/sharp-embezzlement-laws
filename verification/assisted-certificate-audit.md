# Assisted XX certificate audit

This report records an independent AI-assisted mathematical audit accompanying the paper.

The original uniform game separates the central-cut open XX ground state with a shared Φ₉ from every finite-rank maximally entangled input at both tested separating lengths. The exact floor is γ = (3 − 2√2)/243. Providing a flat competitor the same Φ₉ produces a flat state of rank 9D.

| Chain length | Fixed-budget interval for the prescribed strategy |
| ---: | ---: |
| 262,144 | [0.000746210786, 0.000747368515] |
| 524,288 | [0.000701624233, 0.000702783850] |
| 1,048,576 | [0.000661007046, 0.000662385080] |

The first row bounds the loss of the prescribed strategy. The latter two give the claimed upper bounds 0.000703 and 0.000663, below γ > 0.0007060612150362.

- The two physical cross-cut parity blocks have identical singular lists. The sign-corrected even block is the sum of the cosecant and secant terms. Degree-23 dyadic polynomial compression applies to the full finite matrix. Exact Jacobi root brackets, interval weights, and directed Taylor evaluation enclose the weighted quadrature matrix.
- The primary verifier centers each spectral enclosure at the grid point nearest the exact Rayleigh quotient of a proposed vector. It computes the residual in integers and assigns distinct exact eigenvalues via disjoint signed intervals. The Frobenius mass accounts for all unassigned singular values. Sorting after assigning padded model zeros contracts the singular-list distance.
- The normalized model is placed in the physical Schmidt mode bases. The two copies of the one-parity list give squared vector distance \(2[1-\prod_i\cos^2(\theta_i-\widehat\theta_i)]\), including angle zero for omitted model modes. Rearrangement contracts the resulting amplitude distance. Exact integer counts include all degeneracies through the model cutoff. For the unnormalized prefix \(u\), positivity gives the full model defect upper bound \(\|(I-D)u\|^2+2\tau\) with \(\tau=1-\|u\|^2\).
- The defect square root transfers to the physical amplitude list with error at most \(2d\). The ninefold repeated amplitudes on \(\psi_N\otimes\Phi_9\) give exact loss \(2\mathcal E(\psi_N)/243\) for the POVMs prescribed in the physical Schmidt basis. The separate fixed-POVM transfer follows from \(\|K^{1/2}(v-w)\|\le\|v-w\|\) for a game loss effect \(0\le K\le I\).

The original spectral centers allowed solver-dependent drift in the last digits of the distance bound. The final verifier uses exact Rayleigh centering and checks fixed distance budgets before scoring. Full runs with one, two, and four OpenBLAS threads satisfied every published error cap and produced identical certified loss intervals, matching the archive. The appendix exponents match the verifier proposals.

The portable independent verifier uses checked Jacobi eigenpairs, a separate fixed-point kernel approximation, and exact rational products for the Schmidt ranks. It passed with one, two, and four OpenBLAS threads. It bounds the physical canonical losses below 0.000702728 and 0.000662269, and the same-POVM transferred losses below 0.000702858 and 0.000662419. The checked source is `reproducible/strategy/assisted_independent_check.py`.
