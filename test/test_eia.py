import pytest
import responses
import json
from unittest.mock import patch, mock_open
import numpy as np
import os

from src.eia import EIA

# File open mocking modeled on https://stackoverflow.com/questions/1289894/how-do-i-mock-an-open-used-in-a-with-statement-using-the-mock-framework-in-pyth
@pytest.fixture
@patch(
    "builtins.open",
    new_callable=mock_open,
    read_data="BA1,BA2\r\nFirst.Name.1,Second.Name.2",
)
def eia(mock_file, monkeypatch):
    monkeypatch.setenv("EIA_API_KEY", "fakepas")
    e = EIA()
    assert e.key == "fakepas"
    assert e.regions == {"BA1": "First.Name.1", "BA2": "Second.Name.2"}
    mock_file.assert_called_with(e.cache + e.REGION_F, "r")
    return e


@pytest.fixture
def mock_response() -> dict:
    with open("test/external/eia930.json", "r") as f:
        response = f.read()
    return json.loads(response)


"""
Test caching: get First.Name.1, then get it again.
"""


@responses.activate
def test_get_series(eia, mock_response):
    responses.add(
        responses.GET,
        eia.BASE_URL_SERIES.format("fakepas", "First.Name.1"),
        json=mock_response,
    )
    series = eia.get_net_generation("BA1", 2022)
    assert len(series) == 10
    assert series["20220303T05Z"] == 10420

    # Kill endpoint and ask again -- we should read from file.
    responses.add(
        responses.GET,
        eia.BASE_URL_SERIES.format("", "First.Name.1"),
        json="Don't read me!",
    )
    series_read = eia.get_net_generation("BA1", 2022)
    assert np.all(series == series_read)

    os.remove(eia.cache + eia.regions["BA1"])
