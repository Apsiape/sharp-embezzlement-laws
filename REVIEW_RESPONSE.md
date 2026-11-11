# Treatment of the review feedback

Revision 4 integrates the substantive proof, scope and attribution suggestions in the 11 November 2026 feedback. The feedback summary was checked against the manuscript and primary sources; its mathematical verdict was not treated as a substitute for the additional audit recorded in this package.

| Feedback | Treatment |
| --- | --- |
| Unsupported XXZ gapped/split sentence | Removed. The spectral argument does not require it. |
| CTM spectrum imported rather than derived | The abstract, main theorem, corollary, caption and discussion explicitly state the CTM input. The physical identification is separated from the spectral theorem. |
| Magic-square bridge sketched | A lemma states the hypotheses and gives support compression and the matrix-unit transport, with exact citations to the preceding game paper. |
| `2^84` rounding analysis and thin printed margin | The inequality chain is expanded. The public bound is rounded to `0.000461`; the sharper certified endpoint remains recorded. |
| Smoothing sketch loose | The main text points to the atom-inclusive counting proof. The sharp tower remainder is attributed to the direct counting criterion, rather than to the weaker smooth-spectrum transfer. |
| Earlier Gaussian crossover | Zhang–Liu–Chen is credited at the crossover statement. The root-fidelity convention, LOCC setting, and this paper's local-unitary and fixed-target conclusions are distinguished. |
| Earlier maximally entangled floor | Mančinska's flagged-square result is credited, together with earlier Bell-separation work. The explicit constant and certified physical-chain strategies are identified as the additions here. |
| Four-spin tilted-CHSH result | Acín–Massar–Pironio and Bamps–Pironio are credited explicitly. The chain ground state realizes a known separation; the original-game branch is played classically. |
| Infinite-chain and QFT attribution | Matsui, Keyl–Matsui–Schlingemann–Werner, and van Luijk–Stottmeister–Werner–Wilming are credited at the relevant claims. The wedge locator is Theorem 97 of arXiv:2401.07299v4, dated 2025. |
| XXZ entropy and counting input | Calabrese–Cardy, Alba–Calabrese–Tonni, and Calabrese–Cardy–Peschel receive precise references. The symmetric zero-mode and symmetry-broken conventions are distinguished. |
| Earlier rates | Credits added, but the proposed claim that the sharp law “supersedes both” is not used: the earlier works control different quantities, and Kera's stated density estimate is not itself a critical-chain upper rate. |
| Known entropy and variance coefficients | Identified as recoveries of known boundary entropy and capacity laws. |
| Bibliography | Corrected Harris, Otte, Chern and Naaijkens–Ogata entries; included Werner in the classification attribution; cited de Boer et al.; removed the unused Basteiro entry. |
| Source and sample housekeeping | Removed the unused `xx-proof.tex`; made the singlet remainder's parameter dependence explicit; specified the Ising sample size and physical bipartitions. |
| Finite-size interpretation | Added pre-asymptotic interpretation for the displayed numerical sizes, and the exact assisted threshold with certified comparisons. |
| Optional numerical optimization of `7/3` | Omitted: an observed numerical ratio is not a uniform certified replacement for the proved bound. |

## Two precision changes to the feedback itself

The statement that the asymptotic remainder dominates at every computable size is too broad. The revision discusses the displayed finite samples and their measured deviations from the leading term.

Kera et al. bound Gaussian conversion error in terms of covariance spectral density. Sierant studies a dyadic rank-window quantity using a smoothed spectral model. Neither statement is the optimal squared local-unitary error proved here, so the paper compares the settings and conventions without claiming an improved exponent for the same quantity.

The new audits are AI-assisted mathematical reviews. The package reports their scope separately from executable arithmetic certificates and does not claim Lean verification or journal peer review.
