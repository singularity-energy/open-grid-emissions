import io
import zipfile

import pytest
import requests

import oge.download_data as download_data
import oge.constants as constants


YEARS_TO_TEST = list(
    reversed(range(constants.earliest_validated_year, constants.latest_validated_year))
)
HOURLY_YEARS_TO_TEST = list(
    reversed(
        range(constants.earliest_hourly_data_year, constants.latest_validated_year)
    )
)


def test_download_pudl_data():
    """Make sure that PUDL data download works."""
    download_data.download_pudl_data(source="aws")
    print("DONE")


def test_download_egrid():
    """Make sure that eGRID download works."""
    download_data.download_egrid_files()
    print("DONE")


def test_download_chalendar_files():
    """Make sure that Chalendar download works."""
    download_data.download_chalendar_files()
    print("DONE")


def test_download_eia_electric_power_annual():
    """Make sure that we can download electric power annual data from EIA."""
    download_data.download_eia_electric_power_annual()
    print("DONE")


def test_download_eia930():
    """Test that EIA-930 data download works for all years."""
    print("Will test the following years:\n", HOURLY_YEARS_TO_TEST)
    for year in HOURLY_YEARS_TO_TEST:
        print(f"Testing EIA-930 download for {year}")
        download_data.download_raw_eia930(years_to_download=[year])
    print("DONE")


def test_download_epa_psdc():
    """Test that the EPA Power Sector Data crosswalk download works."""
    download_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.3/epa_eia_crosswalk.csv"
    )
    print("DONE")


def test_download_raw_eia860():
    """Test that EIA-860 data download works."""
    print("Will test the following years:\n", YEARS_TO_TEST)
    for year in YEARS_TO_TEST:
        print(f"Testing EIA-860 download for {year}")
        download_data.download_raw_eia860(year)


def test_download_raw_eia923():
    """Test that EIA-923 data download works."""
    testable_years = range(2008, 2021)
    print("Will test the following years:\n", testable_years)
    for year in testable_years:
        print(f"Testing EIA-923 download for {year}")
        download_data.download_raw_eia923(year)


def test_download_raw_eia_906_920():
    testable_years = range(2005, 2008)
    for year in testable_years:
        print(f"Testing EIA-906/920 download for {year}")
        download_data.download_raw_eia_906_920(year)


class FakeResponse:
    """Minimal stand-in for `requests.Response`, failing for URLs that EIA does not
    serve."""

    def __init__(self, url: str, available_urls: list[str]):
        self.url = url
        self.ok = url in available_urls

    def raise_for_status(self):
        if not self.ok:
            raise requests.HTTPError(f"503 Server Error for url: {self.url}")

    def iter_content(self, chunk_size: int):
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w") as zf:
            zf.writestr("data.xlsx", "data")
        yield buffer.getvalue()


def patch_requests(monkeypatch, tmp_path, available_urls: list[str]) -> list[str]:
    """Patch `requests.get` and the downloads folder, returning the list of requested
    URLs."""
    requested = []

    def fake_get(url, **_):
        requested.append(url)
        return FakeResponse(url, available_urls)

    monkeypatch.setattr(download_data.requests, "get", fake_get)
    monkeypatch.setattr(
        download_data, "downloads_folder", lambda path: str(tmp_path / path)
    )
    return requested


def test_download_helper_raises_on_failed_request(monkeypatch, tmp_path):
    """A failed request should raise an error without saving the error page."""
    patch_requests(monkeypatch, tmp_path, available_urls=[])
    download_path = tmp_path / "file.zip"
    with pytest.raises(requests.HTTPError):
        download_data.download_helper(
            "https://www.eia.gov/missing.zip",
            str(download_path),
            str(tmp_path / "file"),
            requires_unzip=True,
        )
    assert not download_path.exists()


def test_download_raw_eia860_early_release_from_archive(monkeypatch, tmp_path):
    """Early release EIA-860 data should be downloaded from the archive, keeping the
    ER suffix, once EIA has moved it there."""
    monkeypatch.setattr(download_data, "current_early_release_year", 2025)
    monkeypatch.setattr(download_data, "latest_validated_year", 2024)
    archive_url = (
        "https://www.eia.gov/electricity/data/eia860/archive/xls/eia8602025ER.zip"
    )
    requested = patch_requests(monkeypatch, tmp_path, available_urls=[archive_url])

    download_data.download_raw_eia860(2025)

    assert requested[-1] == archive_url
    assert (tmp_path / "eia860/eia8602025ER/data.xlsx").exists()


def test_download_raw_eia923_early_release_from_archive(monkeypatch, tmp_path):
    """Early release EIA-923 data should be downloaded from the archive, keeping the
    er suffix, once EIA has moved it there."""
    monkeypatch.setattr(download_data, "current_early_release_year", 2025)
    monkeypatch.setattr(download_data, "latest_validated_year", 2024)
    archive_url = (
        "https://www.eia.gov/electricity/data/eia923/archive/xls/f923_2025er.zip"
    )
    requested = patch_requests(monkeypatch, tmp_path, available_urls=[archive_url])

    download_data.download_raw_eia923(2025)

    assert requested[-1] == archive_url
    assert (tmp_path / "eia923/f923_2025er/data.xlsx").exists()


def test_get_pudl_release_uses_requested_build_and_data_store(tmp_path):
    """The release should be read from the marker file of the requested build in the
    requested folder, regardless of the PUDL_BUILD and PUDL_DATA_STORE settings."""
    for build, marker in [
        ("stable", "stable-2026-10-14-1552-90561c620-v2026-10-0"),
        ("nightly", "nightly-2026-10-20-0602-1a2b3c4d5-main"),
    ]:
        (tmp_path / build).mkdir()
        (tmp_path / build / marker).touch()

    assert (
        download_data.get_pudl_release("stable", str(tmp_path / "stable"))
        == "stable build v2026.10.0 (built 2026-10-14)"
    )
    assert (
        download_data.get_pudl_release("nightly", str(tmp_path / "nightly"))
        == "nightly build main (built 2026-10-20)"
    )
    assert (
        download_data.get_pudl_release("nightly", str(tmp_path / "stable"))
        == "nightly build (unknown release)"
    )
