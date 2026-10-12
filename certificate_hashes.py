"""Explicit, versioned text hashes; historical digests remain immutable.

The one legacy mapping is an audited CRLF serialization of the committed
two-cusp ledger, not a general 'try another hash until it passes' rule.
"""
from __future__ import annotations

import hashlib

TEXT_HASH_POLICY = "sha256-utf8-lf/v1"
LEGACY_POLICY = "track-b-two-cusp-crlf-compatibility/v1"
LEGACY_LEDGER = {
    "383557121919fcddfb8a3a6b8592238f612eca6b1f6bba3bc4583d15611281a4":
    "8ddb66397840d9ebddca7ac2d4e8810076ec031c0931ea4c21cc0940c66c6486",
}


def lf_bytes(raw: bytes) -> bytes:
    raw.decode("utf-8", errors="strict")
    out = raw.replace(b"\r\n", b"\n")
    if b"\r" in out:
        raise ValueError("bare CR is outside sha256-utf8-lf/v1")
    return out


def text_hash(raw: bytes) -> str:
    return hashlib.sha256(lf_bytes(raw)).hexdigest()


def ledger_hash_evidence(raw: bytes, expected: str, policy: str | None = None,
                         schema: str | None = None) -> dict:
    actual = hashlib.sha256(raw).hexdigest()
    evidence = {"raw_sha256": actual, "recorded_sha256": expected,
                "policy": policy or "sha256-raw-bytes/v1", "matches": False,
                "historical_digest_rewritten": False}
    if policy == TEXT_HASH_POLICY:
        evidence["canonical_sha256"] = text_hash(raw)
        evidence["matches"] = evidence["canonical_sha256"] == expected
    elif policy is not None:
        evidence["reason"] = "unknown hash policy"
    elif actual == expected:
        evidence["matches"] = True
    elif schema == "track-b-two-cusp-hejhal/v1" and expected in LEGACY_LEDGER:
        normalized = text_hash(raw)
        evidence.update(policy=LEGACY_POLICY, canonical_sha256=normalized,
                        registered_canonical_sha256=LEGACY_LEDGER[expected],
                        matches=normalized == LEGACY_LEDGER[expected])
    return evidence
