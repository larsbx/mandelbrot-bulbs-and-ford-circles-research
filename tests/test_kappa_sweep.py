import json

import pytest

from kappa_sweep import load_checkpoint


def test_checkpoint_discards_only_malformed_trailing_record(tmp_path):
    checkpoint = tmp_path / "checkpoint.jsonl"
    first = {"p": 1, "q": 7}
    second = {"p": 2, "q": 7}
    complete = (json.dumps(first) + "\n" + json.dumps(second) + "\n").encode()
    checkpoint.write_bytes(complete + b'{"p": 3, "q"')

    assert load_checkpoint(checkpoint) == {1: first, 2: second}
    assert checkpoint.read_bytes() == complete


def test_checkpoint_rejects_malformed_nontrailing_record(tmp_path):
    checkpoint = tmp_path / "checkpoint.jsonl"
    checkpoint.write_text('{"p": 1}\nnot-json\n{"p": 2}\n')

    with pytest.raises(json.JSONDecodeError):
        load_checkpoint(checkpoint)
