"""Manual-edit provenance never changes measurements; it labels objects for stratification."""
import json

import numpy as np
import pytest

from organoid_analysis.quantification.mask_provenance import (
    MaskProvenance,
    labels_sha256,
    load_edit_log,
    parse_ids,
    stratification_note,
)


def _labels():
    labels = np.zeros((4, 10, 10), np.uint16)
    labels[1:3, 1:4, 1:4] = 3
    labels[1:3, 6:9, 6:9] = 8
    return labels


def test_hash_depends_on_content_dtype_and_shape():
    labels = _labels()
    assert labels_sha256(labels) == labels_sha256(labels.copy())
    assert labels_sha256(labels) != labels_sha256(labels.astype(np.int32))
    assert labels_sha256(labels) != labels_sha256(labels.reshape(2, 20, 10))
    edited = labels.copy()
    edited[0, 0, 0] = 3
    assert labels_sha256(labels) != labels_sha256(edited)


@pytest.mark.parametrize("value, expected", [("3;8", {3, 8}), ("3, 8", {3, 8}), ([8], {8}), ("", None), (None, None)])
def test_parse_ids(value, expected):
    assert parse_ids(value) == (frozenset(expected) if expected is not None else None)


@pytest.mark.parametrize("value", ["0", "-2", "a", 5])
def test_parse_ids_rejects_non_label_values(value):
    with pytest.raises(ValueError):
        parse_ids(value)


def test_states_by_declaration():
    labels = _labels()
    assert MaskProvenance.pipeline(labels).state(3) == "no"
    assert MaskProvenance.imported(labels, "none").state(8) == "no"
    assert MaskProvenance.imported(labels, "manual").state(8) == "yes"
    assert MaskProvenance.imported(labels, "").state(8) == "unknown"
    assert MaskProvenance.imported(labels, "").declaration == "not_declared"
    partial = MaskProvenance.imported(labels, "Edited", "8")
    assert (partial.state(3), partial.state(8)) == ("no", "yes")
    assert MaskProvenance.imported(labels, "edited").state(3) == "yes"   # no IDs listed: all flagged


def test_invalid_declarations_are_rejected():
    labels = _labels()
    with pytest.raises(ValueError, match="one of"):
        MaskProvenance.imported(labels, "some")
    with pytest.raises(ValueError, match="only accompany"):
        MaskProvenance.imported(labels, "manual", "3")
    with pytest.raises(ValueError, match="not present"):
        MaskProvenance.imported(labels, "edited", "3;99")


def test_summary_and_note():
    provenance = MaskProvenance.imported(_labels(), "edited", [8])
    summary = provenance.summary([3, 8])
    assert (summary["n_edited_no"], summary["n_edited_yes"], summary["n_edited_unknown"]) == (1, 1, 0)
    assert summary["edited_ids"] == [8]
    assert stratification_note(["no", "no"]) is None
    assert "1 object(s) have manually edited" in stratification_note(["no", "yes"])


def test_edit_log_schema(tmp_path):
    path = tmp_path / "log.json"
    path.write_text(json.dumps({"cell": {"declaration": "edited", "edited_ids": [2]}}))
    assert load_edit_log(path)["cell"]["edited_ids"] == [2]
    path.write_text(json.dumps({"tissue": {"declaration": "none"}}))
    with pytest.raises(ValueError, match="keyed by"):
        load_edit_log(path)
    path.write_text(json.dumps({"cell": {"declaration": "none", "who": "me"}}))
    with pytest.raises(ValueError, match="only"):
        load_edit_log(path)
