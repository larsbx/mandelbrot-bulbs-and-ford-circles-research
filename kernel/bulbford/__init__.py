"""bulbford — satellite-bulb geometry of unicritical families vs Ford circles.

Modules: cf (arithmetic of p/q), dynamics (cycles, multipliers, bulb proxies),
index (parabolic holomorphic index ι_{p/q}), taylor (Taylor data of ρ(u), u = q²ε).

The generic exact arithmetic of p/q and of angle doubling (continued fractions,
modular inverses, units, Farey sequences, rotation cycles, wakes, rotation
numbers, doubling orbits) is `rational_dynamics_py`, and the Krawczyk operator,
its preconditioner and the box tests of certify/antipode are `root_isolation_py`
over `closed_interval`; all three are vendored byte-for-byte from
larsbx/finite-math-kernels into vendor/python and pinned in vendored.toml.
`cf`, `wake` and `cycles` keep their signatures as thin adapters over the first.
Importing this package puts a source checkout's vendor/python first on sys.path
when it is not already there, then refuses any copy of the three other than
the pinned one: the checkout's vendor/python, or the one an installed wheel
carries beside this package.
"""
import importlib
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_VENDOR = _HERE.parents[1] / "vendor" / "python"
if _VENDOR.is_dir() and str(_VENDOR) not in sys.path:
    sys.path.insert(0, str(_VENDOR))

# Each is checked as soon as it is imported, before anything it imports can be read from a shadowing copy.
for _name in ("rational_dynamics_py", "closed_interval", "root_isolation_py"):
    _FOUND = Path(importlib.import_module(_name).__file__).resolve().parent
    if _FOUND.parent not in (_VENDOR, _HERE.parent):
        raise ImportError(f"{_name} resolved to {_FOUND}, not the pinned copy in {_VENDOR} "
                          "(vendored.toml); put that directory first on sys.path")
del _name, _FOUND
