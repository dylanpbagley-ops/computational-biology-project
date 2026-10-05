# Frozen protocol hash-definition erratum

The protocol file was preserved byte-for-byte. Its full-file SHA-256 is 9cf1278621e7b520bfd2bbfbac0a35065bbe60fb4f6b7ddd1f467e75add6b268.

The frozen protocol metadata says its protocol_sha256 covers canonical JSON for all fields except protocol_sha256. That wording is inconsistent with the stored value. Recomputing canonical sorted-key compact UTF-8 JSON while excluding only protocol_sha256 yields c6a95fe8cbd20698872932361ce33eea112922805a69b2f6e956cefce124ed02, not the stored 45103a84c6d99f9bea919f8ae5b3511f0fbf515b1e52b5d4f976cdea6b2167a5. Excluding both protocol_sha256 and protocol_hash_definition yields the stored value 45103a84c6d99f9bea919f8ae5b3511f0fbf515b1e52b5d4f976cdea6b2167a5; this covers all substantive protocol fields while omitting both hash metadata fields.

The validator now checks the original full-file SHA, the stored protocol hash, the single-exclusion result, and the two-metadata-field result explicitly. It also checks that this erratum is present. The protocol, thresholds, inputs, data, and numerical outcomes were not changed. This note documents the inconsistency rather than silently rewriting the frozen protocol or its hash-definition text.
