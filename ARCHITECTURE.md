# Repository architecture

This repository adopts the estate repository template `estate-repository-v1`,
whose canonical source is `larsbx/estate-governance`.

The machine-readable source of repository structure and authority is
[`ESTATE.toml`](ESTATE.toml): this repository's estate position (SPEC_estate v0.1)
and its layout. The contract, `estate-repository-template-v2`, and the audit live
only in `larsbx/estate-governance`; nothing from it is vendored here. CI checks
governance out at the commit the `estate-governance` `[[dep]]` pins and runs the
audit from there, and the audit verifies its own sha256 against that pin.
The ordering rule is:

```text
authority -> mathematical/domain concern -> implementation language
```

Python under `kernel/bulbford/` is the canonical executable. The claim register
`RESEARCH_bulb-ford-correction.md` holds claim state; `EXACT_EVIDENCE_BOUNDARY.md`
separates exact conformance vectors from ball/FFT/continuation evidence;
`reference/legacy/` and `experiments/scripts/` are non-authoritative.

`vendor/python/` is the `pinned_external` plane: `rational_dynamics_py` and
`vendoring` from `larsbx/finite-math-kernels`, copied byte-for-byte and pinned
per file in `vendored.toml`, whose digest the `finite-math-kernels` `[[dep]]` of
`ESTATE.toml` carries. The kernel's generic `p/q` and doubling arithmetic calls
it; nothing in `vendor/` is edited here (README, "Vendoring").

The layout is canonical: every plane in `ESTATE.toml` maps exactly its `target`
(root-level files aside) and no migration step is pending; the audit enforces
both. Directory renames alone must not change claim status, acceptance, or
authority.
