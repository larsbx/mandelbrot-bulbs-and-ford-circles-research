"""Where the sweeps read and write: experiments/data, independent of the working directory.

Importing this also puts vendor/python (the vendored finite-math-kernels packages,
pinned in vendored.toml) on sys.path, so a script run with PYTHONPATH=kernel can
import `rational_dynamics_py`.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "experiments" / "data"
VENDOR = ROOT / "vendor" / "python"
if str(VENDOR) not in sys.path:
    sys.path.append(str(VENDOR))
