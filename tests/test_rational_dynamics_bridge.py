"""Replay the shared rational_dynamics R1 contract in the Ford-bulb consumer."""

import json
from fractions import Fraction
from pathlib import Path

from bulbford.cf import cf, convergent_denominators, from_cf, modinv, xstar


def test_shared_rational_dynamics_r1_vectors():
    path = Path(__file__).resolve().parents[1] / "data" / "rational_dynamics_r1_vectors.json"
    payload = json.loads(path.read_text())
    assert payload["schema"] == "rational-dynamics-r1-fixtures/v1"
    assert payload["source_repository"] == "larsbx/finite-math-kernels"
    assert len(payload["source_commit"]) == 40

    for row in payload["vectors"]:
        p, q = row["p"], row["q"]
        inv = modinv(p, q)
        signed = inv - q if 2 * inv > q else inv
        expansion = cf(p, q)

        assert inv == row["inverse"]
        assert signed == row["signed_inverse"]
        assert expansion == tuple(row["continued_fraction"])
        assert from_cf(expansion) == (p, q)
        assert convergent_denominators(expansion) == tuple(row["convergent_denominators"])
        assert xstar(p, q) == Fraction(abs(signed), q)

        # The shared R1 kernel also owns doubling modulo one. This consumer
        # does not need a second implementation, so replay its finite vector
        # directly with integer arithmetic.
        doubled_num = (2 * p) % q
        doubled = Fraction(doubled_num, q)
        expected = Fraction(row["double_num"], row["double_den"])
        assert doubled == expected

        # For 0 < p < q in canonical CF form, the Ford identity is exactly
        # the shared-kernel identity pinned upstream.
        assert row["convergent_denominators"][-2] == abs(signed)
