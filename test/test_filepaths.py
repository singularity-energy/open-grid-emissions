import pytest

from oge import filepaths


@pytest.mark.parametrize(
    "store, expected",
    [
        (None, False),
        ("local", False),
        ("LOCAL", False),
        ("1", False),
        ("s3", True),
        ("S3", True),
        (" s3 ", True),
        ("2", True),
    ],
)
def test_use_pudl_s3(monkeypatch, store, expected):
    """`PUDL_DATA_STORE` should be read the same way regardless of case."""
    if store is None:
        monkeypatch.delenv("PUDL_DATA_STORE", raising=False)
    else:
        monkeypatch.setenv("PUDL_DATA_STORE", store)
    assert filepaths.use_pudl_s3() == expected
    assert filepaths.get_pudl_data_store().startswith("s3://") == expected


def test_use_pudl_s3_raises_on_unknown_value(monkeypatch):
    """A misspelled value should fail clearly instead of silently picking a store."""
    monkeypatch.setenv("PUDL_DATA_STORE", "aws")
    with pytest.raises(ValueError, match="PUDL_DATA_STORE"):
        filepaths.use_pudl_s3()
