"""bulbford — satellite-bulb geometry of unicritical families vs Ford circles.

Modules: cf (arithmetic of p/q), dynamics (cycles, multipliers, bulb proxies),
index (parabolic holomorphic index ι_{p/q}), taylor (Taylor data of ρ(u), u = q²ε).

The generic exact arithmetic of p/q and of angle doubling (continued fractions,
modular inverses, units, Farey sequences, rotation cycles, wakes, rotation
numbers, doubling orbits) is `rational_dynamics_py`, vendored byte-for-byte
from larsbx/finite-math-kernels into vendor/python and pinned in vendored.toml.
`cf`, `wake` and `cycles` keep their signatures as thin adapters over it.
Importing this package puts vendor/python on sys.path when it is not already
there (pytest's pythonpath and CI's PYTHONPATH put it there first).
"""
import sys
from pathlib import Path

_VENDOR = str(Path(__file__).resolve().parents[2] / "vendor" / "python")
if _VENDOR not in sys.path:
    sys.path.append(_VENDOR)
