# Changes in version 1.0.0 (first public release)

Version 1.0.0 is revision 4 with the wording fixes from a refutation-first check of that revision. The check found no error in the mathematics: the magic-square bridge lemma and the explicit `2^84` rounding chain were rederived, and the chain-only certificate was rerun. No theorem, constant or certificate changes.

- The rank-window criterion of Zanoni, Theurer and Gour is now described as deciding universal *LOCC* embezzlement (Properties section). The paper's error is defined through local unitaries, and Sierant's Section S6.3 shows that the two notions differ.
- The numerical appendix states the enumeration cutoff the code uses: products above relative probability `10^-11` of the largest product; a run that would keep more than 12,000,000 entries stops with an error.
- The comparison with Kera, van Luijk, Stottmeister and Wilming now states their Theorem 1 (for an `ε`-dense covariance spectrum, trace-distance error at most `11 ε^(1/4)` times the number of modes, minimized over Gaussian unitaries) and says that its rate and ours concern different quantities.
- The singlet-family proof states why the infinite-chain ground state is unique.
- The disclosure is now a section titled "How this work was done", naming the AI models used and the independent checks, and separating the exact certificates from the floating-point computations.
- The title page carries the version. The PDF is renamed `sharp-embezzlement-laws.pdf`; licence files and release metadata are added; the checksum manifest is regenerated.

# Changes in revision 4

## Claims and attribution

The title and lead XX/Ising theorem are unchanged. The Gaussian crossover is now credited where stated to Zhang–Liu–Chen's earlier LOCC root-fidelity law; this paper proves the local-unitary chain result and its sharper fixed-target expansion. The entropy and information-variance asymptotics are explicitly identified as recoveries of known entropy and capacity laws.

The qualitative flat-state floor is credited to Mančinska's flagged magic-square result. The paper identifies its contribution as the explicit constant, the quantitative comparison, and the certified chain separations. Liang–Vértesi–Brunner is included among the earlier separation results. Acín–Massar–Pironio receives credit for the tilted-CHSH bound, alongside Bamps–Pironio's sum-of-squares proofs.

The XXZ theorem, abstract, introduction, game-loss corollary, figure caption and discussion now state the result for the imported corner-transfer-matrix spectrum. An unnecessary uncited assertion that the corresponding half-chain bipartition is split is removed. The physical Wick-violation result and exact circuit examples retain their separate hypotheses and proofs.

Matsui's XX type result, the Keyl–Matsui–Schlingemann–Werner formulation, van Luijk–Stottmeister–Werner–Wilming's modular obstruction and wedge-vacuum results, and the prior CHSH comparison are explicitly credited. XXZ entropy and counting inputs receive equation-level references. Bibliographic titles, version numbering and publication details are corrected.

## Proof exposition and numerical interpretation

The magic-square bridge is stated as a lemma with its hypotheses, support compression and explicit transport. The sharp tower remainder is tied directly to the atom-inclusive counting criterion; the smoothed Calabrese–Lefèvre comparison is distinguished from that proof.

The chain-only `2^84` rounding analysis is expanded into explicit inequalities. The public upper bound is rounded outward to `0.000461`; the unchanged certificate gives `0.00046079975444316377…`. The assisted threshold `E < (3−2√2)/2` and certified finite-size comparisons are stated explicitly.

The finite samples are identified as pre-asymptotic, the Ising level-ratio sample size is specified, and the singlet-family bounded remainder is allowed to depend on its parameter. Physical half-chain bipartitions are stated explicitly. An unused duplicate proof file and unused bibliography entries are removed.

## Verification

The revision adds a separate mathematical audit of the sharp law, crossover, CTM interpretation, bridge and robustness comparison, plus primary-source checks of the feedback's attribution claims. Existing computational certificates and exact constants are retained. The release package records its reproduction commands, output and checksums. Lean work has not started.
