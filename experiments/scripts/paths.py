"""Where the sweeps read and write: experiments/data, independent of the working directory."""
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
