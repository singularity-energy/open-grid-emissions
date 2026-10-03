"""Uploads the zipped data release files to Zenodo as a new version of the record.

Run as a plain module, after `output_data.zip_data_for_zenodo(year)` has been run for
each year (see `notebooks/manual_data/prepare_uploads.ipynb`):

    python -m oge.zenodo_upload
    python -m oge.zenodo_upload --publish
    python -m oge.zenodo_upload --sandbox --record-id 123456

This always clears out every file inherited from the previous version before
uploading, rather than only replacing files for years that changed, since a version
bump can revise data for any past year, not just the current one.

Requires a Zenodo personal access token (scopes: `deposit:write`, `deposit:actions`),
created at https://zenodo.org/account/settings/applications/ (or, for `--sandbox`,
https://sandbox.zenodo.org/account/settings/applications/) and saved with no other
content to a gitignored file at the top of the repo:

    .zenodo_api_key            (production)
    .zenodo_sandbox_api_key    (sandbox)

By default the script leaves the new version as an unpublished draft so you can
review it on zenodo.org before it becomes permanent (publishing a deposition cannot
be undone). Pass --publish to publish immediately instead.
"""

import argparse
import datetime
import os
from importlib.metadata import version
from pathlib import Path

import requests

from oge.filepaths import data_folder
from oge.logging_util import configure_root_logger, get_logger

logger = get_logger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]

PRODUCTION_RECORD_ID = "10909484"

MAX_ATTEMPTS = 3
REQUEST_TIMEOUT_SECONDS = 300


def get_base_url(sandbox: bool) -> str:
    """Returns the Zenodo API base URL for the given environment."""
    return "https://sandbox.zenodo.org/api" if sandbox else "https://zenodo.org/api"


def get_token(sandbox: bool) -> str:
    """Reads the Zenodo API token from its gitignored file at the repo root."""
    filename = ".zenodo_sandbox_api_key" if sandbox else ".zenodo_api_key"
    token_path = REPO_ROOT / filename
    if not token_path.is_file():
        raise FileNotFoundError(
            f"Could not find a Zenodo API token at {token_path}. Create a personal "
            "access token (scopes: deposit:write, deposit:actions) and save it, "
            f"with no other content, to {token_path}."
        )
    token = token_path.read_text().strip()
    if not token:
        raise ValueError(f"{token_path} is empty.")
    return token


def request_with_retry(
    session: requests.Session, method: str, url: str, **kwargs
) -> requests.Response:
    """Issues a request, retrying on transient failures with a short backoff.

    Uploads of multi-gigabyte zip files are the main reason this exists: a dropped
    connection partway through a large PUT shouldn't require restarting the whole
    release process.
    """
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            response = session.request(
                method, url, timeout=REQUEST_TIMEOUT_SECONDS, **kwargs
            )
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            if attempt == MAX_ATTEMPTS:
                raise
            logger.warning(
                f"{method} {url} failed ({e}); retrying ({attempt}/{MAX_ATTEMPTS})"
            )


def get_latest_deposition_id(base_url: str, record_id: str, session: requests.Session) -> str:
    """Resolves a record id (any version) to the deposition id of its latest version."""
    response = request_with_retry(
        session, "GET", f"{base_url}/records/{record_id}/versions/latest"
    )
    return str(response.json()["id"])


def create_new_version(base_url: str, deposition_id: str, session: requests.Session) -> dict:
    """Creates a new draft version of the deposition and returns its representation."""
    response = request_with_retry(
        session, "POST", f"{base_url}/deposit/depositions/{deposition_id}/actions/newversion"
    )
    draft_url = response.json()["links"]["latest_draft"]
    return request_with_retry(session, "GET", draft_url).json()


def clear_inherited_files(base_url: str, draft: dict, session: requests.Session) -> None:
    """Deletes every file inherited from the previous version of the draft."""
    for file in draft["files"]:
        logger.info(f"removing inherited file {file['filename']}")
        request_with_retry(
            session,
            "DELETE",
            f"{base_url}/deposit/depositions/{draft['id']}/files/{file['id']}",
        )


def upload_files(bucket_url: str, session: requests.Session) -> None:
    """Uploads every zip file in the local zenodo data folder to the draft's bucket."""
    zenodo_folder = data_folder("zenodo")
    filenames = sorted(f for f in os.listdir(zenodo_folder) if f.endswith(".zip"))
    if not filenames:
        raise FileNotFoundError(
            f"No zip files found in {zenodo_folder}. Run "
            "output_data.zip_data_for_zenodo(year) for each year first."
        )
    for filename in filenames:
        path = os.path.join(zenodo_folder, filename)
        size_mb = os.path.getsize(path) / 1e6
        logger.info(f"uploading {filename} ({size_mb:,.0f} MB)")
        for attempt in range(1, MAX_ATTEMPTS + 1):
            try:
                with open(path, "rb") as f:
                    response = session.put(
                        f"{bucket_url}/{filename}",
                        data=f,
                        timeout=REQUEST_TIMEOUT_SECONDS,
                    )
                response.raise_for_status()
                break
            except requests.RequestException as e:
                if attempt == MAX_ATTEMPTS:
                    raise
                logger.warning(
                    f"upload of {filename} failed ({e}); "
                    f"retrying ({attempt}/{MAX_ATTEMPTS})"
                )


def update_metadata(base_url: str, draft: dict, session: requests.Session) -> None:
    """Bumps the draft's version and publication date to match this release."""
    metadata = draft["metadata"]
    metadata["version"] = version("oge")
    metadata["publication_date"] = datetime.date.today().isoformat()
    request_with_retry(
        session,
        "PUT",
        f"{base_url}/deposit/depositions/{draft['id']}",
        json={"metadata": metadata},
    )


def publish(base_url: str, deposition_id: str, session: requests.Session) -> dict:
    """Publishes the draft. This cannot be undone."""
    response = request_with_retry(
        session, "POST", f"{base_url}/deposit/depositions/{deposition_id}/actions/publish"
    )
    return response.json()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sandbox",
        action="store_true",
        help="Upload to sandbox.zenodo.org instead of production, for a dry run.",
    )
    parser.add_argument(
        "--publish",
        action="store_true",
        help="Publish the new version immediately. Cannot be undone. Default is to "
        "leave it as a draft for manual review.",
    )
    parser.add_argument(
        "--record-id",
        default=None,
        help="Zenodo record id to create a new version of (any version's id works). "
        f"Defaults to the production record ({PRODUCTION_RECORD_ID}). Required when "
        "--sandbox is passed, since sandbox records have no relation to production ones.",
    )
    args = parser.parse_args()

    if args.sandbox and args.record_id is None:
        parser.error("--record-id is required when using --sandbox")
    record_id = args.record_id or PRODUCTION_RECORD_ID

    configure_root_logger()

    base_url = get_base_url(args.sandbox)
    token = get_token(args.sandbox)
    session = requests.Session()
    session.headers.update({"Authorization": f"Bearer {token}"})

    latest_id = get_latest_deposition_id(base_url, record_id, session)
    draft = create_new_version(base_url, latest_id, session)
    logger.info(f"created draft version {draft['id']}")

    clear_inherited_files(base_url, draft, session)
    upload_files(draft["links"]["bucket"], session)
    update_metadata(base_url, draft, session)

    logger.info(f"draft ready for review: {draft['links']['html']}")

    if args.publish:
        result = publish(base_url, draft["id"], session)
        logger.info(f"published: {result['links']['record_html']}")
    else:
        logger.info("skipping publish (pass --publish to publish automatically)")


if __name__ == "__main__":
    main()
