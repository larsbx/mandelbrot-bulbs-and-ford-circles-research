"""Where the sweeps read and write: experiments/data, independent of the working directory.

Importing this also imports `bulbford`, which puts vendor/python (the vendored
finite-math-kernels packages, pinned in vendored.toml) first on sys.path and refuses
any other `rational_dynamics_py`, so a script run with PYTHONPATH=kernel imports the
pinned copy.
"""
from pathlib import Path

import bulbford  # noqa: F401  (puts the pinned vendor/python first)

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "experiments" / "data"
