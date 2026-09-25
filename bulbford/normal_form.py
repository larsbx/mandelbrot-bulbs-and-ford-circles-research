"""Equivariant normal form of f(z) = λ₀z + z², λ₀ = e^{2πip/q}:  f∘H = H∘N,  N(z) = λ₀z(1 + b₁z^q + b₂z^{2q}),
H(z) = z + Σ_{k≥2} H_k z^k with H_{q+1} = H_{2q+1} = 0.  Homological equation at degree k:
  H_k (λ₀ − λ₀^k) = [k>q+1]·(k−q) λ₀^{k−q−1} β₁ H_{k−q} − (H²)_k,   β₁ = (H²)_{q+1},  β₂ = (H²)_{2q+1},  b_i = β_i/λ₀.
Composing N q times: a = [z^{q+1}] f^q = q b₁,  ι_{p/q} = (q²−1)/(2q) + b₂/(q b₁²).  O(q²) operations."""
from __future__ import annotations
from typing import NamedTuple
import mpmath as mp


class NormalForm(NamedTuple):
    p: int
    q: int
    H: tuple            # H_0..H_{2q+1}
    b1: mp.mpc
    b2: mp.mpc

    @property
    def a(self) -> mp.mpc:
        return self.q * self.b1

    @property
    def iota(self) -> mp.mpc:
        q = self.q
        return mp.mpf(q * q - 1) / (2 * q) + self.b2 / (q * self.b1 ** 2)


def normal_form(p: int, q: int, dps: int | None = None) -> NormalForm:
    with mp.workdps(dps or q + 40):
        lam = mp.expjpi(mp.mpf(2 * p) / q)
        H = [mp.mpc(0)] * (2 * q + 2); H[1] = mp.mpc(1)
        conv = lambda k: mp.fsum(H[i] * H[k - i] for i in range(1, k))
        beta1 = mp.mpc(0)
        for k in range(2, 2 * q + 1):
            if k == q + 1:
                beta1 = conv(k); continue
            rhs = -conv(k) + ((k - q) * lam ** (k - q - 1) * beta1 * H[k - q] if k > q + 1 else 0)
            H[k] = rhs / (lam - lam ** k)
        beta2 = conv(2 * q + 1)
        return NormalForm(p, q, tuple(H), beta1 / lam, beta2 / lam)
