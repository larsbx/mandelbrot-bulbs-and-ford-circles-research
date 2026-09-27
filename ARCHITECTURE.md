# Repository architecture

This repository adopts the estate repository template `estate-repository-v1`,
whose canonical source is `larsbx/estate-governance`.

The machine-readable source of repository structure and authority is
[`estate.toml`](estate.toml). The reusable contract is
[`docs/architecture/estate-repository-template-v1.md`](docs/architecture/estate-repository-template-v1.md),
vendored byte-for-byte and pinned by sha256 in the `[governance]` table of
`estate.toml`. Do not edit it or `tools/audit_estate_layout.py` locally; change
them in governance and re-vendor.

The ordering rule is:

```text
authority -> mathematical/domain concern -> implementation language
```

Python under `kernel/bulbford/` is the canonical executable. The claim register
`RESEARCH_bulb-ford-correction.md` holds claim state; `EXACT_EVIDENCE_BOUNDARY.md`
separates exact conformance vectors from ball/FFT/continuation evidence;
`reference/legacy/` and `experiments/scripts/` are non-authoritative.

The layout is canonical: every plane in `estate.toml` maps exactly its `target`
(root-level files aside) and no migration step is pending; the audit enforces
both. Directory renames alone must not change claim status, acceptance, or
authority.
