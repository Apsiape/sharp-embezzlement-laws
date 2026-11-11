#!/bin/sh
set -eu
cd "$(dirname "$0")"
unset PYTHONOPTIMIZE
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -B verify_manifest.py
cd reproducible
python -B test_reproducibility.py
python -B certificate.py
python -B constant_verify.py
python -B strategy/check_rounding_bounds.py
python -B witness.py
python -B witness.py --compact
python -B strategy/check_four_site_hybrid.py
python -B strategy/assisted_xx_certificate.py
python -B strategy/assisted_independent_check.py
python -B strategy/certify_xx84_modes.py
python -B strategy/model_score.py strategy/model84_exponents.json --cutoff 380000 --distance-bound 3/2500 --require-separation
