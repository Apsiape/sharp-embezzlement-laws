Programs and data for "Reading the von Neumann type of matter with a finite Bell game" (draft, 10 November 2026).

data/   JSON files with the numbers behind Figures 1(b) and 2 and Tables 2 and 3
        (S in nats, b = -ln p_1, E = embezzlement error, varX = Var(log2 T)).
code/   Python programs (numpy, scipy, matplotlib).
        dilation_defect.py, lib.py, numlib.py   spectra of free-fermion chains and the functional E
        run_crit.py, xx_from_ising.py           critical XX and Ising chains
        run_gapped.py                           XX chain with a staggered potential
        cl_and_pairs.py, towers2.py, towers3.py theory curves and copies of two-qubit states
        validate.py, check_props.py             checks against exact diagonalization and of the inequalities
        make_figs2.py                           figures and table rows

The programs were run in one working directory and contain the absolute path of that directory
(sys.path.insert lines and the path in lib.py). Replace these paths before running them elsewhere.
Note: the function H_functional returns 2*E; E is obtained by dividing by two (done in numlib.stats).
