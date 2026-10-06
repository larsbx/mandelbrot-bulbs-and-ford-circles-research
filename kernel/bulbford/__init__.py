"""bulbford — satellite-bulb geometry of unicritical families vs Ford circles.

Modules: cf (arithmetic of p/q), dynamics (cycles, multipliers, bulb proxies),
index (parabolic holomorphic index ι_{p/q}), taylor (Taylor data of ρ(u), u = q²ε).

The generic exact arithmetic of p/q and of angle doubling (continued fractions,
modular inverses, units, Farey sequences, rotation cycles, wakes, rotation
numbers, doubling orbits) is `rational_dynamics_py`, vendored byte-for-byte
from larsbx/finite-math-kernels into vendor/python and pinned in vendored.toml.
`cf`, `wake` and `cycles` keep their signatures as thin adapters over it.
Importing this package puts a source checkout's vendor/python first on sys.path
when it is not already there, then refuses any `rational_dynamics_py` other than
the pinned copy: the checkout's vendor/python, or the one an installed wheel
carries beside this package.
"""
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_VENDOR = _HERE.parents[1] / "vendor" / "python"
if _VENDOR.is_dir() and str(_VENDOR) not in sys.path:
    sys.path.insert(0, str(_VENDOR))

import rational_dynamics_py as _rd  # noqa: E402

_FOUND = Path(_rd.__file__).resolve().parent
if _FOUND.parent not in (_VENDOR, _HERE.parent):
    raise ImportError(f"rational_dynamics_py resolved to {_FOUND}, not the pinned copy in {_VENDOR} "
                      "(vendored.toml); put that directory first on sys.path")
del _rd, _FOUND
