import hashlib
import unittest
from pathlib import Path
from certificate_hashes import (TEXT_HASH_POLICY, LEGACY_POLICY,
                               ledger_hash_evidence, text_hash)


class HashPolicyTests(unittest.TestCase):
    def test_new_policy_is_portable_and_content_sensitive(self):
        self.assertEqual(text_hash(b'{"x":1}\n'), text_hash(b'{"x":1}\r\n'))
        self.assertNotEqual(text_hash(b'{"x":1}\n'), text_hash(b'{"x":2}\n'))
        self.assertNotEqual(text_hash(b'{"x":1}\n'), text_hash(b'{"x":1} \n'))
        with self.assertRaises(ValueError):
            text_hash(b'a\rb')

    def test_registered_history_is_explicit_and_immutable(self):
        raw = Path('track_b_hejhal_rows.jsonl').read_bytes()
        expected = '383557121919fcddfb8a3a6b8592238f612eca6b1f6bba3bc4583d15611281a4'
        result = ledger_hash_evidence(raw, expected, schema='track-b-two-cusp-hejhal/v1')
        self.assertTrue(result['matches'])
        self.assertEqual(result['policy'], LEGACY_POLICY)
        self.assertFalse(result['historical_digest_rewritten'])
        self.assertFalse(ledger_hash_evidence(raw + b' ', expected,
                         schema='track-b-two-cusp-hejhal/v1')['matches'])
        self.assertFalse(ledger_hash_evidence(raw, expected,
                         schema='unregistered-schema')['matches'])

    def test_unknown_policy_or_unregistered_hash_fails(self):
        raw = b'a\n'
        expected = hashlib.sha256(b'a\r\n').hexdigest()
        self.assertFalse(ledger_hash_evidence(raw, expected)['matches'])
        self.assertFalse(ledger_hash_evidence(raw, text_hash(raw), 'unknown')['matches'])
        self.assertTrue(ledger_hash_evidence(raw, text_hash(raw), TEXT_HASH_POLICY)['matches'])
