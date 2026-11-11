# Independent audit of the finite XX certificate

The $N=2^{84}$ chain-only certificate proves the claimed strict separation. The exact XX covariance, quadrature error, binary64 residual bounds, Gaussian-state comparison, and integer score give a physical loss below $0.000460800$. The number $0.000410720762$ bounds the score on the geometric model before the vector-distance transfer. The audit found no change to the claimed constants necessary.

## Reproduction

From the workspace root, all three commands exited zero with Python 3.12.14, NumPy 2.3.5, and SciPy 1.17.0.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python output/revision3/reproducible/strategy/certify_xx84_modes.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python output/revision3/reproducible/strategy/model_score.py output/revision3/reproducible/strategy/model84_exponents.json --cutoff 380000 --distance-bound 3/2500
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/revision3/audit_physical_certificate.py
```

| Quantity | Verified value or upper bound |
| --- | ---: |
| Rounded covariance dimension | 1418 |
| Retained modes per parity | 81 |
| Omitted squared singular mass | $3.456356338642579\times10^{-11}$ |
| Retained squared singular error | $4.723883784770301\times10^{-11}$ |
| One-parity singular-list error used | $4.81\times10^{-5}$ |
| Aligned vector error from the mode verifier | $0.0011063$ |
| Aligned vector error used in scoring | $3/2500=0.0012$ |
| Schmidt prefix rank | 185,926,736 |
| Tail probability upper endpoint | 0.000082217162715639803344159365716926293998049 |
| Prefix duplication error upper endpoint | 0.013735796367426328639093403676098135191003108 |
| Prefix nine-block variance upper endpoint | 0.000800858721842588365885236909330938907160582 |
| Model loss upper endpoint | 0.000410720761968614898973037738961275050608514 |
| Physical loss upper endpoint | 0.000460799754443163771129846265247324951980458 |
| Tracial floor lower endpoint | 0.000706061215036254742372932311031291534405992 |

The table of 81 increasing exponents agrees with the JSON file. The audit script verifies the two reported mode-error thresholds on the exact `Fraction` returns of the production residual verifier. Its independent score uses one copy of $1+z^k$ per rank update, exact signed polynomials for mixed-block variance, and repeated outward multiplication on a $2^{-200}$ grid to enclose powers of $q$. Its separate upper endpoints for the model and physical losses are 0.000410720761968614898973037738961275050309702 and 0.000460799754443163771129846265247324951663952.

## Exact spin-chain covariance

Put $N=4m$, $n=2m$, $M=N+1$. The open hopping matrix has sine eigenvectors $\sqrt{2/M}\sin(\pi kx/M)$ with eigenvalues proportional to $\cos(\pi k/M)$. There is no zero level. Filling $k=1,\ldots,n$ or its complement gives the same cross-cut singular values. The finite cosine sum is

$$
 \sum_{k=1}^n\cos(\pi k w/M)=\frac{\sin(\pi w/2)}{2\sin(\pi w/(2M))}-\frac12.
$$

For sites $x=n+1-u$, $y=n+v$, the product-to-sum arguments are $y-x=u+v-1$ and $x+y=M+v-u$. The constant terms cancel. Opposite parities of $u,v$ give zero in both numerators. At $u=2i,v=2j$, the first sine numerator is $-(-1)^{i+j}$ and the second is $(-1)^{i+j}$. The minus between them gives the same sign on both terms, and the second denominator is a cosine. A row and column sign change gives the symmetric positive-entry kernel in the appendix.

The alternating-site unitary sends the full covariance $C$ to $I-C$. Its left restriction yields $C_A-I/2=\left(\begin{smallmatrix}0&D\\D^*&0\end{smallmatrix}\right)$ in parity order. Since $C^2=C$, the cross block times its adjoint is $C_A(1-C_A)=I/4-(C_A-I/2)^2$. The two diagonal blocks have equal spectra since $D$ is square and $DD^*,D^*D$ have the same eigenvalues, including multiplicities of zero. This proves the parity duplication.

## Quadrature and compression

For a bin of length $L$, set $F_j(x)=x^{\underline j}(x-L)^{\underline j}$. The degree-$j$ monic polynomial is $j!\Delta^jF_j/(2j)!$. The zeros of $F_j$ on $0,\ldots,j-1$ and $L,\ldots,L+j-1$ remove boundary terms in summation by parts. For the norm, $\nabla^j$ of the monic polynomial equals $j!$, and $(-1)^jF_j(x)=(j!)^2\binom{x}{j}\binom{L-x+j-1}{j}$. Vandermonde's identity gives $h_j=(j!)^4\binom{L+j}{2j+1}/(2j)!$. The squared off-diagonal in the centered coordinate scaled by $1/L$ is

$$
 \frac{h_j}{L^2h_{j-1}}=\frac{j^2(L^2-j^2)}{4L^2(4j^2-1)}.
$$

The full Jacobi matrix is the multiplication operator on the bin. The $(1,1)$ entry of its $p$-th power agrees with that of its $k\times k$ compression for $p<2k$, since a path takes at least $2k$ steps to leave the first $k$ indices and return. The quadrature weights are the squared first eigenvector coordinates times $L$. The nodes are within the bin by the Rayleigh bounds for compressed multiplication. When $k=L$ the quadrature is the original measure.

The degree-17 kernel approximation is a polynomial of separate degree at most 17 in each index on every bin rectangle. Evaluation on the original grid and weighted evaluation on quadrature nodes have the same inner products for those polynomials, by moments through degree 34 or by the exact original measure on short bins. The two approximate matrices act through isometric copies of the same coefficient spaces and share their nonzero singular values. The 1418 quadrature nodes consist of 31 from bins of lengths 1 through 16, 18 each from the 77 larger nonsingleton bins, and one from the final singleton.

If the two bins start at $a,b$ with lengths $L_a,L_b$, the midpoint of $w=2i+2j-1$ is $w_0=2(a+b)+L_a+L_b-3$ and its deviation is at most $L_a+L_b-2$. As $L_a\le a,L_b\le b$, the ratio is at most $1/3$. The exact degree-17 geometric remainder is bounded by $3^{-18}/(\pi w)$. Summing by $i+j-1$ bounds the principal Frobenius norm by $\sqrt{1+\ln(2m)}/(2\pi)<4/3$. The quadrature norm is at most twice the original grid norm by the binwise ratio of extreme $w$ and conservation of bin weight.

For either regular kernel term, write $t_0$ for the real center and $\delta$ for the maximal deviation on its bin rectangle. Both argument intervals are symmetric, and $\delta\le1/4$ because the largest nonsingleton bin has length $m/2$. If $R_0=1$ for $f$ and $R_0=1/2$ for $g$, then $\delta\le R_0-|t_0|$. For the Cauchy disk of radius $R-|t_0|$ about $t_0$, the Taylor ratio is at most $(1/4)/(R-R_0+1/4)$. This is $1/3$ for $R=3/2$ and $1/2$ for $R=3/4$. On the first circle the sine modulus is at least $1/\sqrt2$, by either the imaginary-part hyperbolic sine or the real-part sine, and $|f|\le 1/\sqrt2+2/(3\pi)<1$. On the second circle the cosine modulus is at least $\cos(3\pi/8)>3/8$, hence $|g|<4/3$. Cauchy estimates give the asserted remainders and the sum of the two Frobenius errors is bounded by

$$
 4\cdot3^{-18}+\tfrac12\big[(3/2)3^{-18}+(8/3)2^{-18}\big]
 < 5.099\times10^{-6}<6\times10^{-6}.
$$

Here $1/2$ dominates $2m/M$ because each measure has mass $m$. The singular-value perturbation inequality on both sides of the shared polynomial spectrum gives the claimed list error.

## Rounding and singular-value matching

- **Jacobi.** The rational square comparisons bound each rounded off-diagonal entry within $2\times10^{-15}$, so the symmetric tridiagonal perturbation has norm at most $4\times10^{-15}$. The checked squared norm and residual imply exact normalized residual below $4\times10^{-14}$. Each candidate interval contains an eigenvalue, and the $k$ disjoint intervals exhaust the $k$-dimensional spectrum. The exact off-target gap from a candidate is at least $10^{-3}-4\times10^{-14}$. After sign alignment, the normalized eigenvector distance is at most $\sqrt2$ times the residual divided by that gap. Including normalization gives distance below $5.66\times10^{-11}$. The checked first-coordinate square is above $10^{-3}$, which gives relative weight error below $2\times10^{-7}$ after the two multiplications. Each exact node is at least the bin start $a$, both $a,L$ are binary64-exact, and $L\le a$. The eigenvalue error contributes less than $4\times10^{-14}a$, and forming the midpoint and node adds less than $5\cdot2^{-53}a$. The relative node error is below $5\times10^{-14}$.
- **Kernel.** For the sum argument, node errors give relative error at most $(4/3)5\times10^{-14}$ before arithmetic since $2i+2j-1\ge3$. For the difference argument their absolute contribution is at most $(\pi/2)5\times10^{-14}$ before arithmetic. The denominator $2M$ rounds with relative error $1/M$ and the Machin enclosure bounds the absolute error in the rounded $\pi$ by $10^{-15}$. Including the remaining operations gives the claimed $3\times10^{-13}$ bounds. The Taylor remainders are below $10^{-17}$ by the first omitted terms at 1.58 for sinc and 0.8 for cosine. The absolute coefficient sums at these bounds are below 3. The Horner term paths each have at most 31 rounding factors, counting the argument square and coefficient rounding, so $3\gamma_{31}+10^{-17}<10^{-13}$. The exact factors satisfy $\sin(t)/t\ge2/\pi$ and $\cos(t)\ge1/\sqrt2$ on their respective intervals. The Lipschitz bounds and the polynomial errors give relative errors below $10^{-12}$ in the evaluated sine and cosine. The rounded square and division in the weight-root check add two unit-roundoff factors, so the exact squared ratio to the rounded product differs from one by less than $1.1\times10^{-14}$. Weight relative error, the square-root check, and the remaining positive operations give entrywise relative error below $10^{-6}$. The rounded Frobenius norm is below $\sqrt{4+3\times10^{-12}}$, so the Frobenius difference is at most $10^{-6}\|\widehat K\|_F/(1-10^{-6})<2.1\times10^{-6}$.
- **Trace and eigenpairs.** With unit roundoff $u=2^{-53}$ and $d<1600$, the trace error is at most $4\gamma_{2d+1}/(1-\gamma_{2d+1})<3\times10^{-12}$. The candidate squared-norm test and the positive-sum error imply $\|v\|>.99$. The norm of the arithmetic error in each matrix-vector product is at most $\gamma_{d+1}\|\widehat K\|_F\|v\|$. The computed residual norm is at most $40\cdot10^{-14}$. The scalar multiplication and final subtraction add another unit-roundoff factor each. The exact normalized residual is below $3\times10^{-12}$ with ample slack. The additive rounding bound $2^{-1075}$ per scalar operation for gradual underflow leaves the bounds valid. The signed candidate gaps exceed $10^{-10}$, so their exact eigenvalue intervals of radius $3\times10^{-12}$ are disjoint.
- **Mass and pairing.** The squared Frobenius norm is the sum of squared singular values. Subtracting the squared lower absolute endpoints of the selected distinct exact eigenvalues bounds the mass of all the others. On $0\le r\le1$, $r/(1+r^2)$ is increasing, so the directed Decimal power endpoints give valid rational endpoints for the model singular values. Assign the selected values in the proposed pairing and all unselected values to padded zeros. Sorting both lists decreases or preserves squared distance by rearrangement. This works for every selected set of distinct exact eigenvalues. The exact Fraction checks verify both manuscript mass bounds and also that their sum is below $(9.05\times10^{-6})^2$. The proof uses $4\times10^{-5}$ for the rounded matrix. Together with the two approximation errors this gives $4.81\times10^{-5}$ for the sorted physical singular values. The maximum model singular value is below $0.498337276$, and adding the list error gives a physical maximum below $0.499$.

## Physical vector and score

Wick factorization in the local orbitals that diagonalize $C_A$ gives a product reduced density matrix, with one-mode probabilities $\lambda_j,1-\lambda_j$. The left Jordan–Wigner operators are supported in the left spin block, so these are the eigenvalues of the spin reduced state. Choose the more probable mode occupation as the zero label and use the exact physical Schmidt bases for the model. The one-mode amplitude angle is $\theta(s)=\tfrac12\arcsin(2s)$ and $\theta'(s)=(1-4s^2)^{-1/2}<16$ on $[0,.499]$. The frozen physical modes enter the comparison with model angles zero. For the parity-doubled lists the overlap is the product of squared cosines. Using $1-\prod_i(1-x_i)\le\sum_i x_i$ and $\sin^2 y\le y^2$ bounds the vector distance by $16\sqrt2$ times the one-parity list error. The code uses $23(4.81\times10^{-5})=0.0011063<3/2500$.

The rank recurrence, tail, duplication, and block variance refer to unnormalised prefix amplitudes $q^e\sqrt{p_0}$. Each mode pair contributes exact coefficients $1,2,1$ at exponents $0,k,2k$. Updating from high to low reads the lower indices from the previous factor. Runs record the exact number of Schmidt positions of each exponent. The duplication overlap sweeps the direct runs and the doubled runs, with $1/\sqrt2$ from $D$. Since $D$ is an isometry, the two squared norms in $E$ contribute twice the prefix mass. The variance of a block is its sum of squares minus the square of its sum divided by nine. Constant complete blocks contribute zero. A partial final block uses zero padding. Every interval primitive encloses its exact operation by integer floor or ceiling, including the square root via `isqrt` on scaled integer endpoints.

The local Schmidt-prefix projections yield a common branch on the aligned model. The homogeneous operational bound applies to the prefix, and the classical branch contributes $\tau/27$. The measurements extend locally on the exact physical spaces with any fixed off-support outcomes, using the chain as the only shared resource. For the same measurements on either vector, the losing effect $F$ is a positive contraction since it averages sums of subsets of complete product POVMs. Thus $\|F^{1/2}\psi_N\|\le\|F^{1/2}\widehat\psi\|+\|\psi_N-\widehat\psi\|$. The physical loss bound squares this inequality after inserting $3/2500$, and the exact score places it below $0.000460800$ and the tracial floor.

## Exposition

Four steps can be made explicit in the appendix. State the bin-center Taylor ratios using the symmetric argument intervals. Give the rearrangement argument that passes from the verified pairing to the sorted singular lists. Identify the local Schmidt-prefix projection and classical complement when applying the operational bound. Give the positive-contraction transfer on the same POVMs.
