# Independent audit of the explicit robustness constant

This report records an independent AI-assisted mathematical audit accompanying the paper.

**Verdict: valid, with no proof repair required.**

The final manuscript proves

\[
 \mathcal E(\psi)\le \frac{2\,706\,291}{2}
       \bigl(1-\omega_G(\psi)\bigr)
 =1\,353\,145.5\bigl(1-\omega_G(\psi)\bigr)
\]

for every unit vector in a tensor-product bipartite Hilbert space, arbitrary local POVMs with the game's finite outcome sets, and finite or countably infinite Schmidt rank. The proof does not depend on local dimension, a lower bound on a sector's norm, attainment of the game value, or a finite-rank approximation of the strategy.

This was a fresh audit of the final manuscript sources. I did not read the earlier research audits, review ledgers, or prior verdicts. I checked the operator arguments before checking their numerical certificates. No manuscript file was edited.

## Sources and scope

The audited proof is manuscript/sections/robust-proof.tex, lines 1–370. I also checked the game rules, the definition and Schmidt formula in sections/setting.tex, the comparison statement in sections/properties.tex, and the main theorem in type-of-matter.tex. The upper operational comparison and the other physical applications are outside this verdict.

The source hashes at completion were:

| Source | SHA-256 |
|---|---|
| sections/robust-proof.tex | dfc888dd8715caff4dbd308dbe9073c411ac459a2fd0625761633c4781d9c358 |
| sections/game-rules.tex | dc3297d8e680b47b5a0e5e7a07b018d788906d292a7d088ba29b37fc4c3165ae |
| sections/setting.tex | 481eb6f9de7758bac9edadcc6e08628ae57871991b0f0d83bb23f6618a10c336 |
| sections/properties.tex | 8db347d5cbff03441534cef48bce3ac0c5eb53afa04af354ca8097dc4fc9716c |
| ../reproducible/constant_verify.py | 50ca7d4cbf668967334afd9016eac36154a589e29aea011e0e2d93a76197a3aa |

## 1. Game probabilities and the common projective dilation

There are exactly \(9\cdot25+18=243\) equally weighted supported question pairs. In one copy the internal forbidden mass is

\[
 N=M+T+S:
 \quad 9\text{ row-column pairs},\quad
 4\text{ calibration pairs},\quad
 4\text{ sign pairs},
\]

with the remaining eight internal pairs free. Thus these quantities are sums of conditional losses, not losses already divided by 25 or 243.

The projective reduction in lines 295–307 is a single common local dilation for all questions. For each question the POVM isometry \(W:\mathcal H\to\mathcal H\otimes\mathbb C^m\) satisfies \(W^*W=I\), and

\[
 U_W=\begin{pmatrix}0&W^*\\ W&I-WW^*\end{pmatrix}
\]

is self-adjoint and obeys \(U_W^2=I\). Conjugating an outcome-coordinate PVM by \(U_W\) gives compression \(M_a\) on the same input copy of \(\mathcal H\). Padding to a common finite \(m\) is legitimate. The input inclusion is independent of the question. Applying the two local inclusions preserves every correlation and all nonzero Schmidt amplitudes, including at infinite rank.

After dilation, a binary flag mismatch is exactly
\(\|(F_i-G_j)\psi\|^2\). It is bounded by the full conditional loss \(u_{ij}^2\). Similarly, assigning every flag sign \(+1\) gives
\(\|(A_{ij}-B_{ij})\psi\|\le2u_{ij}\). A one-flag event can contribute to this sign mismatch, but it is already forbidden; this does not invalidate the bound.

## 2. Magic-square parity identities and word rewriting

The flag correction has the correct sign:

\[
 B_{33}=-B_{13}B_{23}+2G_3.
\]

On the nonflag subspace Bob's third-column product is \(-I\); on the flag subspace all three observables are \(+I\). There is no missing factor in this correction.

I checked all three displayed anticommutator expressions on each party:

\[
 d\{b,c\}a,\qquad b\{a,d\}c,\qquad c\{a,d\}b.
\]

Each is within \(W=\sum_{ij}\beta_{ij}\), with \(\beta_{ij}=2u_{ij}\), of \(2G_3\psi\). Each of the nine cell errors is used once. For example, on Alice,

\[
 dbca\psi \rightsquigarrow
 a'c'b'd'\psi=B_{31}B_{32}\psi
 \rightsquigarrow A_{33}\psi
 \rightsquigarrow B_{33}\psi,
\]

using the four generator bounds and the three bounds for cells 31, 32, 33. The other term satisfies \(dcba=cdab=A_{23}A_{13}\) exactly and transfers to \(B_{13}B_{23}\) with the remaining two bounds. The other five cases have the same allocation; the exact local commutations are \(ab=ba,cd=dc\) on Alice and \(a'c'=c'a',b'd'=d'b'\) on Bob. Reversing the pair order uses an actual row or column commutation, not an unproved cross-context commutation.

The reversed-word transfer rule is valid for arbitrary words of local unitaries. For a difference \(L\) of two signed words, \(\|L\|\le2\), so a one-letter suffix costs at most \(2\beta_s\). This is the coefficient used in the rewriting table.

The flag errors in that table are also valid without assuming that Bob's different column flags commute:

\[
 \|(I-G_1)G_3\psi\|
 \le \|(I-G_1)(G_3-F_1)\psi\|
       +\|(I-G_1)F_1\psi\|
 \le u_{13}+d_A,
\]

\[
 \|(I-F_1)G_3\psi\|\le u_{13}.
\]

For an anticommutator without a suffix, the stated trailing-letter transfer commutes with the opposite anchor flag: use \(a'\) or \(c'\) against \(G_1\) on Alice's side, and \(a\) or \(b\) against \(F_1\) on Bob's side. For the \((a,d)\) swap already carrying a suffix, the displayed identity directly supplies the necessary suffix. This verifies all twelve nonzero cases of the table, including the cases where an unjustified projected suffix transfer could otherwise lose a term.

The group-register maps are isometries because each evaluated word is unitary. Their registers have central sign \(-1\). The exact intertwinings for Alice's \(a,b\) and Bob's \(a',c'\) follow from the final two commuting letters in their chosen normal orders. For the other generators the four presence/absence cases of the last two letters have equal register weight \(1/4\), and the empty case is exact. Restricting to
\(\mu=(I-F_1)(I-G_1)\psi\) adds \(2d_A\) or \(2d_B\), as stated.

## 3. Consistency errors and the extraction coefficient

The four consistency estimates on \(\mu\) have the stated coefficients:

\[
\begin{aligned}
 \delta_a&\le2u_{11},&
 \delta_b&\le2u_{12}+2d_B,\\
 \delta_c&\le2u_{21}+2d_A,&
 \delta_d&\le2u_{22}+u_{12}+u_{21}+d_A+d_B.
\end{aligned}
\]

For \(d\), the two comparisons with \((I-G_2)\psi\) and \((I-F_2)\psi\) each cost twice their distance. Taking their arithmetic mean gives the displayed bound. No anchoring commutation for \(d\) with \(F_1\) or \(G_1\) is assumed.

In the central \(-1\) sector, the representation is the real two-qubit Pauli representation tensored with multiplicity. The joint \(a,b\) eigenspaces are related by the commuting \(d,c\) flips, giving this decomposition directly. The four two-party involutions \(K_h\) commute, and their common \(+1\) subspace is exactly \(\Phi_4\) tensored with the two local multiplicity spaces. Thus

\[
 4\eta^2\le\sum_h\|(I-K_h)w\|^2
 \le\sum_h(\alpha_{Ah}+\delta_h+\alpha_{Bh})^2.
\]

The weights \((1,6),(1,5),(2,6,5)\), together with the \(1/2\) root-mean-square factors, give exactly the manuscript's coefficients

\[
 1,\;7,\;6,\;13/2,\qquad
 7/24,\;3/10,\;13/24,\;13/20.
\]

I rebuilt the four coefficient matrices directly from the rewriting table, using explicit integer rows, rather than importing them from the verifier. In coordinate order
\((u_{11},u_{12},u_{13},u_{21},u_{22},u_{23},u_{31},u_{32},u_{33},d_A,d_B)\),
the matrices are:

    L_Ac:
    (4,0,0,4,0,0,0,0,0,2,0)
    (6,2,4,2,2,2,2,2,2,4,0)
    (10,6,4,6,2,2,2,2,2,4,0)

    L_Ad:
    (2,2,4,6,2,2,2,2,2,4,0)
    (0,4,0,0,4,0,0,0,0,2,0)
    (2,6,4,2,6,2,2,2,2,4,0)

    L_Bb:
    (4,4,0,0,0,0,0,0,0,0,2)
    (6,2,4,2,2,2,2,2,2,0,2)
    (10,6,4,6,2,2,2,2,2,0,2)

    L_Bd:
    (2,6,4,2,2,2,2,2,2,0,2)
    (0,0,0,4,4,0,0,0,0,0,2)
    (2,2,4,6,6,2,2,2,2,0,2)

Exact rational arithmetic reproduces every displayed certificate ceiling:

\[
 (521543,521579,521308,521552,521430,
 521308,521308,521308,521308,521546,521180).
\]

The largest actual comparison ratio is
\(9748309/18690<521.6\).
Because \(H\) is entrywise nonnegative, the entrywise inequality with the positive vector \(q\) does imply \(H\le521.6D_*\) as quadratic forms, by the stated pairwise weighted inequality. Also

\[
 v^TD_*v=.561u_{11}^2+\sum_{ij\ne11}u_{ij}^2
       +.439(d_A^2+d_B^2)\le M.
\]

Consequently \(\eta^2\le130.4M\), and \(\eta\le11.42\sqrt M\) is valid.

## 4. Comparison of the unnormalized sectors

No sector is normalized or divided by its norm anywhere in this argument. This is essential when the flag or magic sector is small.

The flag/extraction embeddings are local isometries. The three error components have mutually orthogonal supports and norms \(d_A,d_B,\eta\). Exact intertwining identifies Alice's \(R_1\) outcomes and Bob's \(C_1\) outcomes with the stated \(Z\)- and \(X\)-basis projections.

The sign support expansion in lines 225–227 is an exact operator identity. The two calibration errors contribute at most \(2\sqrt T\) in the four-component direct sum. The third term contributes \(2\sqrt{d^2+\eta^2}\). For the final term, the flag block is counted four times and

\[
 \sum_{x,y}R_yC_xR_y=|0\rangle\langle0|\otimes I\le I
\]

on the magic register. It therefore contributes
\(\sqrt{\eta^2+4d_A^2+d_B^2}\). The apparent asymmetry between \(d_A\) and \(d_B\) is correct.

The selected magic-register overlap is

\[
 (C_x\otimes R_y)\Phi_4
 =\frac{(-1)^{xy}}{2\sqrt2}|c_xr_y\rangle.
\]

Thus both the sign \(t_{xy}=(-1)^{xy}\) and the scale \(h=j/(2\sqrt2)\) are correct. The four sign losses give \(4S\) exactly, and multiplication of the sign-support error by the operator-norm bound 2 gives

\[
 B\le2\sqrt S+4\sqrt T+
 4\sqrt{d^2+\eta^2}+2\sqrt{\eta^2+4d_A^2+d_B^2}.
\]

I checked the four-vector comparison rather than assuming it as a self-testing lemma. Projection onto \(PQ\) gives the common-space CHSH differences for the flag compressions. Projection onto \(C_xR_y\), followed by removal of the fixed qubit factors and \(t_{xy}\), gives the same four-vector pattern for the magic compressions on \(h\). The cycle identity bounds each required compression norm by the sum of these four squared differences with coefficient one.

The four orthogonal components of \(U_0\chi_{00}-V_0\varphi\) are:

1. the \(PQ\) component, controlled by the flag cycle;
2. the \(C_0^A\) component, controlled by the magic cycle;
3. \(P(I-Q)b_{00}\);
4. \((I-P-C_0^A)R_0^B b_{00}\).

The last two are orthogonal to both cycle projections of \(b_{00}\). Consequently no squared error is counted twice, and
\(\|U_0\chi_{00}-V_0\varphi\|\le B\), with no extra \(\sqrt2\) factor.
Schmidt contraction and local unitarity then give
\(\|\operatorname{Sch}(h)-\operatorname{Sch}(\varphi)\|_2\le B\).

The fourfold list of \(w_0=\Phi_4\otimes j\) is
\(4\times(\sqrt2\,\operatorname{Sch}(h))\), whose list scaling has norm multiplier \(2\sqrt2\). Including the two embedding errors gives precisely

\[
 \sqrt{d_B^2+\eta^2}+2\sqrt2d_A+2\sqrt2B
 \le209.6\sqrt M+8\sqrt2\sqrt T+4\sqrt2\sqrt S.
\]

The rational upper bound used for the \(M\) coefficient is
\(13098650677/62500000=209.578410832<209.6\).
Finally,
\(210^2-(209.6^2+160)=196/25>0\).
This proves the \(210\sqrt N\) flagged-sector estimate, including every operator-to-scalar factor.

## 5. Polar label gluing and Schmidt-list alignment

The \(F_j\) across different copies need not commute or be orthogonal. The polar argument does not assume either property.

For \(T v=\sum_jF_jv\otimes|j\rangle\), one has \(A=T^*T=\sum_jF_j\). Extending the polar partial isometry on \(\ker A\) into one unused copy of the original space gives an isometry \(V\) satisfying \(T=V\sqrt A\). The extension works in infinite dimensions as well.

Since \(A\ge F_j\), square-root monotonicity gives \(\sqrt A\ge F_j\). The orthogonal Bob labels imply the blockwise inequality
\(\sqrt A\ge\sum_jF_jQ_j\). Thus

\[
 \|(V-T)\psi\|^2
 \le\langle A+I-2\sum_jF_jQ_j\rangle
 =\sum_j\|(F_j-Q_j)\psi\|^2=N_L.
\]

The raw gluing error after Bob's label isometry has orthogonal register components, and its squared norm is exactly the manuscript's sum of \(\|Q_j(F_k-Q_k)\psi\|^2\). It is at most \(N_L\). Restricting the input by any Bob label projection cannot increase either error, because the Alice maps commute with that projection. The total error is therefore \(2\sqrt{N_L}\) for every subset, with the same common isometry.

In the target the sectors are orthogonal on both local sides, so its Schmidt list is the sorted union of their lists. Replacing a flag by its label costs \(d_j\) for the complementary sector and \(2\sqrt2d_j\) for the fourfold rescaled flag list. This gives exactly

\[
 \left\|\bigsqcup_{k\ne j}S_k-4\times(\sqrt2 S_j)\right\|_2
 \le210\sqrt{N_j}+(1+2\sqrt2)d_j+2\sqrt{N_L}.
\]

The direct sum of the nine comparisons and the \(\ell^2\) triangle inequality give \(210\sqrt{N_F}+(7+2\sqrt2)\sqrt{N_L}\). Sorting is contractive for nonnegative \(\ell^2\) lists. Every \(S_k\) occurs eight times on the left and four times on the right, giving \(2\Delta(U)\), not \(\Delta(U)\) or \(4\Delta(U)\).

The all-label comparison gives \(\|\Lambda-U\|_2\le2\sqrt{N_L}\). Duplication has norm \(\sqrt2\), so \(\Delta\) is \(2\sqrt2\)-Lipschitz. Therefore

\[
 \Delta(\Lambda)\le105\sqrt{N_F}
      +\left(\frac72+5\sqrt2\right)\sqrt{N_L}.
\]

The exact squared label coefficient is
\(249/4+35\sqrt2<112\).

## 6. Final constant and infinite-rank scope

The budgets \(N_F\) and \(N_L\) use disjoint question pairs. The second family of label tests is not needed in this estimate; omitting its nonnegative loss is legitimate:

\[
 N_F+N_L\le243\varepsilon.
\]

The Schmidt formula gives \(\mathcal E=\Delta(\Lambda)^2/2\). Weighted Cauchy–Schwarz now yields

\[
 \mathcal E(\psi)
 \le\frac{105^2+112}{2}(N_F+N_L)
 \le\frac{243(210^2+448)}8\varepsilon
 =\frac{2\,706\,291}{2}\varepsilon.
\]

Taking the infimum over strategies is valid even when the game value is not attained.

The infinite-rank extension has no unresolved limit:

- All physical and extraction operators are bounded, all group and question sums are finite, and all estimates are homogeneous in the input vector.
- Every vector has a countable square-summable Schmidt list. The Schmidt contraction follows from the Hilbert–Schmidt singular-value inequality, with finite-rank approximation if desired.
- Sorted union and decreasing rearrangement are well-defined for these finite unions of square-summable lists. Their contraction extends from finite lists by \(\ell^2\) approximation.
- The polar extension has an explicitly available unused local summand, so it does not require a finite-dimensional rank argument.
- The Schmidt formula for \(\mathcal E\) needs only a supremal local-unitary overlap in infinite dimensions. Alignment of finite Schmidt prefixes approaches the required overlap; exact unitary matching of infinite supports is unnecessary.
- None of these extensions changes any coefficient or normalizes a vanishing sector.

## Checks performed and issues

I ran the supplied constant_verify.py successfully and independently reconstructed its Gram certificate from explicit coefficient rows with rational arithmetic. The arithmetic check supports the constants; the operator derivations above are the proof audit.

**Substantive issues found: none. Required repairs: none.**

The most plausible failure points—third-column parity, projected suffix transfers, the four-cycle orthogonality accounting, nonorthogonal Alice label flags, and infinite-rank local-unitary alignment—are resolved by the arguments actually present in this manuscript. The result should retain the numerical constant \(1\,353\,145.5\), without adding a finite-rank assumption or a sector-normalization assumption.
