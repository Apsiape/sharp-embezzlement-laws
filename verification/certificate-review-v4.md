# Revision 4 finite-chain certificate check

## Chain-only certificate

The appendix gives the binary64 error inequalities for the existing
1418-dimensional covariance calculation. The Jacobi residuals certify distinct
roots and their aligned vectors. With the two multiplications in the weight
calculation and the explicit node arithmetic, the relative weight and node
errors are below $2\times10^{-7}$ and $5\times10^{-14}$. The kernel proof
bounds both argument errors, the Horner truncation and rounding errors, and the
verified square root of each weight product. It uses positivity of the two
reciprocal terms to conclude entrywise relative error below $10^{-6}$.

The squared Frobenius-norm error is below $3\times10^{-12}$ and the exact
normalized eigenpair residuals are below $3\times10^{-12}$, including an
additive allowance for gradual underflow. The resulting Frobenius perturbation
bound $2.1\times10^{-6}$ combines with the quadrature and singular-matching
bounds to give one-parity error $4.81\times10^{-5}$. On the verified singular
range $[0,0.499]$ the Schmidt-angle derivative is below 16. The two equal parity
lists give aligned physical vector distance at most
$16\sqrt2\,(4.81\times10^{-5})<0.0011063<3/2500$.

The two prefix projections give a common branch in the aligned model. The
classical complement contributes its mass divided by 27. The same local
measurements define a physical losing effect $0\le F\le I$, so the triangle
inequality bounds the physical square-root loss by the model square-root loss
plus $3/2500$. The archived model score upper endpoint is
`0.000410720761968614898973037738961275050608514`. An exact integer square root
and outward squaring on the verifier's $10^{-45}$ grid reproduce the physical
loss upper endpoint `0.000460799754443163771129846265247324951980458`. The
manuscript reports the rounded theorem bound $0.000461$.

The revision 3 independent audit in `verification/physical-certificate-audit-v3.md`
checks the exact covariance, compression and physical transfer and independently
recomputes the model score. The two chain-only verifier scripts and exponent
list are byte-identical to the audited revision. The editing audit reproduces
the transfer from the archived upper endpoint. The portable check transfers
the rounded model cap $0.000410720762$ to the stated physical cap $0.000461$.

## Assisted certificate

The prescribed strategy on $\psi_N\otimes\Phi_9$ has loss
$2\mathcal E(\psi_N)/243$. It beats the tracial floor exactly when
$\mathcal E(\psi_N)<(3-2\sqrt2)/2=0.0857864376\ldots$.

| XX sites | Rounded numerical estimate for $\mathcal E(\psi_N)$ | Outward certified interval |
| --- | ---: | --- |
| $2^{18}$ | 0.0908 | [0.0906646, 0.0908053] |
| $2^{19}$ | 0.0854 | [0.0852473, 0.0853883] |
| $2^{20}$ | 0.0805 | [0.0803123, 0.0804798] |

The certified intervals resolve the threshold in all three rows. The
four-decimal estimates in the main text come from the separately labeled
numerical estimates in the archived output. The printed intervals come from
the output's certified physical defect field. The loss intervals and 24 model
exponents in each appendix row match the unchanged verifier and reference JSON.

The four-spin section credits the tilted CHSH construction and the sharp
score to Acín, Massar and Pironio and cites Bamps and Pironio for the
sum-of-squares decompositions. The assisted theorem refers to the explicitly
attributed operational floor.

## Commands and results

Run from the extracted release root.

```sh
python reproducible/strategy/check_rounding_bounds.py
```

The portable check passed the binary64 and compression inequalities, verified
the rounded chain-only loss transfer, and reconstructed the assisted score
intervals and threshold comparisons. It also passed with `python -I -S` in a
folder containing the script and the packaged assisted JSON.

During editing, an additional check reproduced the sharper archived
chain-only endpoint and verified the manuscript model exponents and outward
intervals against the reference data. The staged certificate revision compiled
without warnings, undefined references, or overfull boxes. Extracted PDF text
confirmed the threshold, full chain-only endpoint and assisted intervals.
The final release build and complete executable run are recorded in
`VERIFICATION.md` and `verification-output.txt`.
