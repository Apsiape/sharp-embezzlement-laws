# Independent mathematical audit of the v4 revision

Date: 11 November 2026. No manuscript source was edited by this auditor. No Lean work was started.

## Verdict and scope

**No theorem-breaking error was found in the audited sharp finite-chain law, its stated second-order remainder, the uniform growing-target crossover, or the exact tower law.** The arguments operate on the exact discrete Schmidt lists, including their multiplicities. The stronger remainder does not follow merely from a Gaussian approximation; the microscopic count and second-difference arguments actually supply it.

The v3 XXZ statement needed a scope correction. The final v4 text makes the theorem and its game corollary statements about the explicitly defined bipartite state with the imported CTM spectrum. That correction is sufficient for the spectral proof audited here. It does **not** establish the CTM identification as an operator-algebra theorem for the boundary-selected XXZ states, and the revised text no longer purports to do so.

I additionally checked the standalone proof of the constant C = 2,706,291/2, including its operator reductions, rewriting table, sector comparison, polar gluing, and an independent rational reconstruction of its Gram certificate. I found no repair needed. I checked the new v4 linked modular bridge against the specified game source and verified the standard-pair extraction algebra and modular signs.

The constant-C verdict is for pure states in tensor-product bipartite Hilbert spaces with arbitrary finite-outcome local POVMs, including finite or countably infinite Schmidt rank. It is not an extension of that Schmidt-list argument to general nonfactorizing commuting algebras.

| Claim | Independent finding |
|---|---|
| Finite even open XX, middle cut: E_d = (log d)^2/(4V) + O_d((log N)^-2) | Proof checks out from the exact covariance and discrete counting criterion. |
| Finite even critical open Ising, middle cut | Exact Majorana doubling transfers the XX estimates with the stated factor of two. |
| Uniform growing-target crossover | The uniformity over all integer d is justified by strong L2 convergence of square roots and a uniform translation-operator bound. |
| Exact odd towers, fixed r and d | The O(epsilon^2) conclusion follows directly from the counting criterion, not from the weaker smooth-density transfer alone. |
| XXZ CTM spectrum | Correct as a theorem about the supplied spectrum and canonical state psi_t. The physical identification remains an input. |
| Standalone robustness constant C | Operator proof and rational arithmetic checked; no missing sector-normalization or finite-rank assumption found. |
| New linked modular bridge | Correct under the stated normality/faithfulness hypotheses and cited single-copy bridge; exact extraction is correctly restricted to the standard pair. |

This is **not** a certification of the whole paper. I did not independently audit the arbitrary-cut extension, physical numerical certificates, optimal flat-spectrum game floor, operational upper strategy in full, classical-value enumeration, singlet construction, ADW Wick-violation proposition, interacting circuit/impurity theorems, or the complete infinite-chain/QFT spectral arguments. I did not independently settle priority or all bibliographic attributions. The external review's global verdict is not the basis for this report.

## Material inspected

The stable initial target was `output/revision3/manuscript/type-of-matter.tex` and these sections:

- `law.tex`, `xx-full-proof.tex`, `xx-mesh.tex`, `rank-variance.tex`, `tower-proof.tex`;
- the CTM parts of `interacting.tex` and `xxz-proof.tex`;
- `robust-proof.tex`, with `game-rules.tex` and the Schmidt formula in `setting.tex`.

The finite-law proof was examined before reading the external feedback's conclusion. For the constant C, I also read `research/revision3/constant_cold_audit.md`; the checks below rederive the operator steps and arithmetic rather than treating that earlier verdict as evidence.

After the editors reported completion I checked the current v4 `law.tex` diff, the CTM theorem/proof/corollary/caption, relevant abstract/introduction/discussion passages, and `setting.tex` lines 67–118. For the last check I read Lemmas 1.6, 2.1–2.2 and Section 2.4 of `research/revision4/sources/infinite/type-games-v1.1.0.pdf`.

## 1. The physical finite-chain input

### Hilbert-to-Carleman comparison

The ordered l1 mesh in `xx-mesh.tex` is a sufficiently strong input; it is not just a density-of-states asymptotic.

1. The difference between H_a and H_2 is sign definite because its quadratic form is an integral against the fixed-sign weight t^(a-1) - t. Its trace norm is therefore the displayed convergent sum of diagonal differences, uniformly in the section size.
2. The finite Carleman operator is positive trace class. Its unit-cell compression dominates H_2, since the extra Laplace-integrand factor is at least one. Compression/min–max gives a coordinatewise eigenvalue comparison, so the l1 spectral difference is exactly the bounded trace difference. No unjustified trace-norm comparison of a compression with the whole operator is used here.
3. Under x = e^s/2 the kernel becomes (2 cosh((s-t)/2))^-1. Periodization has the stated Fourier eigenvalues pi sech(2 pi^2 k/L).
4. The potentially troublesome image adjacent to an endpoint becomes the Hankel kernel K(s+t). Pairing consecutive alternating exponential terms makes its rank-two trace norms O(j^-2). More distant images have exponentially summable trace norms. This really does give a uniform trace-norm bound for periodization; an operator-norm or Hilbert–Schmidt bound would not have sufficed.

### Exact covariance and active energies

The sine-covariance cross block splits into the displayed principal Hilbert blocks after parity separation. The analytic remainders admit rank-one Taylor expansions with summable nuclear norms. The two block shifts are 1/2 and 3/2, and their scale is 1/(2 pi); these yield the two mesh branches with singular values (1/2) sech(2 pi^2 k/L_i).

For a projection covariance, B B* = C_A - C_A^2. Hence the cross-block singular values are sqrt(lambda(1-lambda)), with their full multiplicities, and the one-mode energy is 2 arcosh(1/(2s)). On a fixed positive energy band this inverse map is Lipschitz. Bounded l1 singular-value discrepancy therefore removes only O(1) band entries and leaves total matched energy discrepancy O(1).

The tilted Bernoulli variance is bounded below on that band. For frequencies away from zero, its mesh sine-square sum has lower bound cW - C(1+|t|). The O(1) total energy discrepancy costs only another O(|t|), so a sufficiently small fixed tau gives exponential damping up to |t| = tau W. This step remains valid even though the exact finite energies are not an arithmetic progression.

For entropy and information variance, the summands are Lipschitz functions of s all the way to 0 and 1/2. The two branches have total positive-energy density (L_1+L_2)/(2 pi^2), and each relevant integral equals pi^2/6. The resulting S and V are (log N)/6 + O(1), with the stated constants. The separate Holder estimate for the largest-weight summand gives b_N = (log N)/12 + O((log N)^(1/3)), sufficient for the displayed strategy corollary's p_1 = o(E_2) premise.

The critical Ising transfer is an exact finite Majorana decomposition: the XX chain on 2L sites is two uniform Majorana paths, and when L is even the middle cut restricts each path to its first L Majoranas. The Ising energy list is therefore repeated twice in the XX list. This checks the factor 1/12 in S and V and the factor 3/log L in E_d. Odd L and general cuts were not audited here.

## 2. The counting criterion and the O(W^-2) remainder

The proposition's analytic/tilted hypotheses are substantive. Leading entropy or Renyi asymptotics alone would not imply its conclusion. The physical proof above supplies the additional hypotheses rather than assuming them.

### Atoms and strict cumulative counts

In `xx-full-proof.tex` lines 15–75, the sinc^4 smoothing kernel is nonnegative, has finite second moment, and has the claimed compact Fourier support. The assumed damping gives an O(W^-1/2) bound on the smoothed density. Positivity of the kernel near zero then gives

    P_beta(A in I) <= C (|I| + Delta) / sqrt(W),
    Delta = 1/(tau W),

for closed intervals and for singletons. Thus an individual atom has tilted mass O(W^-3/2).

The counting kernel exp(-beta u) 1_{u >= 0} has total variation two. Integrating its shift difference against the preceding concentration estimate bounds the smoothing error by O(Delta/sqrt(W)). This argument explicitly includes its jump at zero.

In the central Fourier integral the cubic cumulant term and the linear kernel term are odd. After rescaling by sqrt(v_beta), the remaining integrated error is relative O(W^-1). The absolute smoothing/atom error has exactly this relative size against the leading W^-1/2 term. Consequently both R(E) and R(E-) have logarithm F(E) + O(W^-1), including at occupied energies. This is a genuine microscopic estimate; no passage from a weak CLT to individual atoms is hidden here.

### Inverse ranks and smooth profile

For the exact step density p(x) = p_ceil(x), each plateau satisfies R(E-) < x <= R(E). Both logarithmic count bounds therefore trap the actual plateau energy within O(W^-1) of F^-1(log x). The inverse derivative is uniformly bounded because F' = beta + O(W^-1). This controls every rank inside a degenerate plateau.

The original tilted law has sub-Gaussian concentration on the cutoff window B = sqrt(W) log W. The law of the exact rank energy under g(y)^2 dy is exactly the Gibbs law of A. The ratio estimate on the central window and the superpolynomial tail bound give ||g-G||_2 = O(W^-1), including finite-rank tails.

For the score a(E) = 1 - 1/F'(E), the derivative calculation has the correct sign and orders:

    a(mu) = O(W^-1),
    a'(mu) = -1/V + O(W^-2),
    a'' = O(W^-2), a''' = O(W^-3).

The moment calculation is the important improvement over a crude Cauchy–Schwarz estimate. E[U] = 0 and E[U^3] = O(W) control the constant-linear and quadratic-linear cross terms; E[U^4] = O(W^2) controls the Taylor remainder. No residual log W factor is needed. It gives

    ||G'||_2^2 = 1/(4V) + O(W^-2),
    ||G''||_2 = O(W^-1).

### Why the stronger remainder survives the approximation

The rank-coordinate identity is exact for integer d:

    E_d = ||(I-T_h)g||_2^2,  h = log d.

It holds because p(x/d) = p_ceil(j/d) almost everywhere on (j-1,j]. Normalization is exact as well, including the interval 0 < x <= 1.

Writing r = g-G, the proof correctly uses

    ||D_h g||^2 - ||D_h G||^2
      = 2 Re <D_h* D_h G, r> + ||D_h r||^2.

The second-difference multiplier is bounded by h^2 t^2. Thus the cross term is O(h^2 W^-2), not the weaker O(W^-3/2) obtained by bounding D_h G and D_h r separately. The fourth-order cosine remainder is bounded by a constant times h^4 t^4 for all real h,t. Together these give precisely

    |E_d - h^2/(4V)| <= C (1+h^2+h^4)/W^2,

with C independent of d. This proves the fixed-d sharp remainders stated for XX and Ising. I found no circular import: the microscopic-counting and profile lemmas use the proposition's hypotheses, not its conclusion.

## 3. The uniform growing-target crossover

The proof produces pointwise convergence of the rescaled **exact** rank density to the standard normal density. Both are normalized. Dominated convergence for their minimum yields L1 convergence, hence L2 convergence of their square roots.

For every real shift a, ||I-T_a|| <= 2. Thus the difference between the squared translation errors of the two square roots tends to zero uniformly over **all** shifts, not only bounded a. Taking a = (log d)/sqrt(V) proves the claimed uniform-in-d absolute o(1) crossover, including d growing far beyond the crossover scale.

The Gaussian square-root overlap is exp(-a^2/8). Therefore E_d tends to zero if and only if log d = o(sqrt(V)). In that regime the polynomial quantitative bound above supplies the relative asymptotic E_d ~ (log d)^2/(4V): after dividing by the leading term, the error tends to zero since log d >= log 2 and (log d)^2 = o(W).

This proof is stronger than a weak convergence assertion. There is no unjustified claim that weak Gaussian convergence alone controls fixed microscopic rank shifts.

The independent `rank-variance.tex` lower bound also checks out. Its density q integrates to one, the exponential tail bound for log(q/p) is valid, log^2(e+y) is concave for y >= 0, and the two shifted means give the displayed translation lower bound. It provides a lower limit only and is not used circularly to establish the matching upper limit.

## 4. Towers and continuum comparisons

The odd-level Fourier bound excludes an unaccounted parity resonance. At theta close to pi the odd energies give a sine-square sum of order M, not a second major arc. The saddle expansion has an odd imaginary cubic correction, and the integrated remainder has the claimed relative O(n^-1/2) size. The cumulative count follows by a controlled sum of the positive coefficient asymptotics.

The uniform inverse-count bound |H^-1(x)-n| <= K_r is justified by comparing the O(n^-1/2) logarithmic count error with the derivative of log H, also of order n^-1/2. It gives the stated uniform logarithmic ratio and O(epsilon) Hellinger comparison to the smooth density, with all degeneracies retained.

However, that smooth comparison together with the smooth-law proposition gives only an O(epsilon^(3/2)) error by the displayed Hellinger transfer estimate. The sharper O(epsilon^2) tower result is correctly proved later, at `tower-proof.tex` lines 392–407, by applying the general counting criterion directly to the discrete physical energy. This is not a logical flaw in the theorem, but the v3 explanatory paragraph described the weaker route as if it proved the sharp remainder. I reported that mismatch to the editor. The final v4 paragraph now correctly names the direct criterion route and treats the Calabrese–Lefevre comparison as additional information.

The smooth Calabrese–Lefevre distribution itself is normalized exactly, including its top interval. Its energy variable is compound Poisson with the displayed moments. The treatment of the derivative jump at y=0 is correct: the proof does not falsely place that profile in H2. These auxiliary results are compatible with, and not substitutes for, the discrete theorem.

## 5. XXZ: what has and has not been established

For the supplied integer-level spectrum, the product and modular identity are correct:

    log Z(t) = pi^2/(12t) - (log 2)/2 + t/24
               - log Z(2 pi^2/t).

Consequently S = pi^2/(6t) - (log 2)/2 + O(t^-1 exp(-2 pi^2/t)) and V = pi^2/(6t) plus an exponentially small correction of the stated order. With t=2 eta, these reproduce all three displayed constants and remainders.

The revised analytic-logarithm estimate is valid in Re z > 1/2. The first floor(1/t) modes give v_beta >= c/t. Their sine-square sum gives the required damping up to pi/t, uniformly in the fixed real tilt interval. Hence the O_d(eta^2) extraction error is a legitimate application of the same criterion.

Adding the zero-energy binary mode duplicates every weight after division by two. On the exact rank axis this translates g by log 2. Translation commutes with D_h and preserves its norm, so the error is exactly unchanged for every integer target dimension, while entropy increases by log 2.

The final v4 theorem defines the normalized canonical bipartite CTM vector psi_t. The strategy corollary is also explicitly on psi_t, and its proof uses an exponentially small p_1 to obtain the stated leading loss. The abstract, introduction, discussion, plot caption, and later circuit comparison now retain that scope. The unsupported v3 assertion that the actual phases are gapped and split is absent from this theorem.

This resolves the scope problem without claiming a new proof of the physical CTM identification. The separate Wick-violation proposition concerns the ADW boundary states and is outside this audit; it should not be included in a blanket endorsement of the CTM theorem.

## 6. Independent check of the explicit robustness constant

I rechecked the proof against the exact game predicates, including the distinction between unnormalized sums of forbidden conditional masses and the final factor 243.

- The common projective dilation is one fixed input inclusion for all questions and preserves Schmidt amplitudes. The displayed block operator is a self-adjoint unitary.
- The third-column parity correction is +2G_3. Each of the three anticommutator identities uses the nine cell errors once. The reversed-word rule and its one-letter suffix cost are correct.
- I checked all twelve nonzero cases of the rewriting table. The projected anticommutator bounds use a transferred letter that actually commutes with the opposite anchor flag. Cases with an existing suffix use the corresponding identity directly. No unproved commutation among different flags is assumed.
- The four consistency bounds on the doubly unflagged vector follow with the stated coefficients. For d_0 the two distance estimates are averaged; there is no missing factor of two.
- The finite-group extraction maps are isometries, the central -1 sector is the two-qubit representation with multiplicity, and the four consistency involutions have joint +1 space Phi_4 times multiplicity. The weighted Cauchy–Schwarz coefficients in H follow exactly.
- In the sign-support expansion, the flag and magic components have the claimed orthogonality. The four-cycle estimate for the two sector vectors counts the leftover pieces only once. The coefficient 2 sqrt(2) in the fourfold Schmidt-list rescaling is correct.
- No flag or magic sector is normalized or divided by its norm. Vanishing or very small sectors therefore cause no hidden assumption.
- The label-gluing argument uses polar decomposition and square-root monotonicity, not orthogonality of Alice's nine flags. The unused summand permits extension to an isometry in infinite dimensions. Sorted-union contraction and the 2 sqrt(2)-Lipschitz estimate for Delta give the stated final label coefficient.

I independently rebuilt the rational Gram calculation from explicit integer rows, without importing the verifier's matrix. Its ceilings were

    521543, 521579, 521308, 521552, 521430, 521308,
    521308, 521308, 521308, 521546, 521180.

The maximum exact ratio was 9748309/18690 < 521.6. The remaining arithmetic gives

    210^2 - (209.6^2 + 160) = 196/25,
    C = (243/2) (105^2 + 112) = 2706291/2.

The proof uses only bounded operators and finite sums before the standard Hilbert–Schmidt Schmidt-list contraction. Its infinite-rank extension needs a supremal unitary overlap, not attainment. I found no circular dependence on the imported robust theorem: that theorem is credited, while this constant is proved directly.

## 7. New v4 linked bridge and exact extraction wording

The revised statement in `setting.tex` is consistent with the source game paper's normal support reduction and calibrated single-copy bridge.

Compression to the closure of M_A Omega is faithful and normal when the Alice marginal is faithful. The subspace reduces M_A, and compressing Bob's effects preserves both POVM normalization and every correlation, even if Bob did not preserve that subspace originally. A perfect binary match makes the Alice effect a centralizer projection. The label tests identify the final magic projection with F_(j+1)+F_(j+2).

Transport by a matrix-unit partial isometry in the centralizer preserves the modular eigenvalue. The sign is correct: sigma_t(w_j) = 2^(it) w_j implies phi(w_j w_j*) = 2 phi(w_j* w_j). The maximum principle gives all nine flag masses 1/9.

For the now-explicit standard pair, s_r = sum_j w_j* F_(j+r) satisfies

    s_r* s_s = delta_rs 1,
    sum_r s_r s_r* = 1,
    sigma_t(s_r) = 2^(-it) s_r.

Tomita's identity gives (J s_s* J) Omega = 2^(-1/2) s_s Omega. Therefore the two column maps have coefficient s_r*(J s_s* J)Omega = delta_rs Omega/sqrt(2), exactly as written. The change removes the old implication that a faithful marginal alone identifies the original Bob algebra with the full standard commutant. The modular-eigenvector obstruction does not require that stronger identification.

This check accepts the explicitly cited normal single-copy bridge as an input after verifying the actual source statement. It is not an independent reproof of every theorem in the game source.

## 8. Source integrity and remaining release constraints

The following load-bearing v4 appendices are byte-for-byte identical to the stable v3 files audited above:

| Section | SHA-256 |
|---|---|
| `xx-full-proof.tex` | `c18e3741d58ba7009be90a74092848747060fe45fe36d47256b2dea8c503e8f4` |
| `xx-mesh.tex` | `2180e834a834208e70a74e4f639801117e7110cb2ba055884a00ac00e23904f8` |
| `rank-variance.tex` | `45b42c231c5b6faa887a874f948ae2bbcb09c95f05a487473305eeb50462ba27` |
| `tower-proof.tex` | `ffe08cef5754a3778b3947a6bfcd67a735944806bb1d26d13da4b925db1d8e4e` |
| `robust-proof.tex` | `dfc888dd8715caff4dbd308dbe9073c411ac459a2fd0625761633c4781d9c358` |

I reviewed the v4 changes in `law.tex`; they preserve the proof and clarify rather than weaken the microscopic hypotheses. I also reviewed the completed CTM and bridge changes as described above.

**Required mathematical repairs remaining in this audited scope: none.** Preserve the CTM input qualification, the exact-discrete versus smooth-profile distinction, and the standard-pair qualification on exact bridge extraction. Any subsequent substantive changes to these proofs require a new check. The unaudited applications and certificates must retain their separate verification records; this report should not be described as a global certification of them.
