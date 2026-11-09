#!/usr/bin/env python3
"""Publish a new version of the CPSR example output record on Zenodo.

The README and pkgdown site link to the concept DOI, which always resolves to
the latest version, so this script can run independently of the version bump.

Standard library only. Three stages, each stopping before the next:

  check    Local validation of the report files. No network access.
  prepare  Run the checks, then create a Zenodo draft with the new files and
           metadata, and verify it. Nothing becomes public.
  publish  Re-verify the draft, ask for typed confirmation, publish, and
           verify the public record. Files of a published version are permanent.

Typical use, on the sandbox first (it has its own account, token and records):

  export ZENODO_SANDBOX_TOKEN=...
  python3 data-raw/zenodo_release.py prepare --sandbox --new-record \\
      --version 2.4.0 --report-dir <dir>
  python3 data-raw/zenodo_release.py publish --sandbox --draft-id <id> \\
      --version 2.4.0 --report-dir <dir>

Then for real, with ZENODO_TOKEN set (scopes: deposit:write, deposit:actions):

  python3 data-raw/zenodo_release.py check   --version 2.4.0 --report-dir <dir>
  python3 data-raw/zenodo_release.py prepare --version 2.4.0 --report-dir <dir>
  python3 data-raw/zenodo_release.py publish --draft-id <id> \\
      --version 2.4.0 --report-dir <dir>
"""

import argparse
import datetime
import gzip
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

CONCEPT_RECORD_ID = "4050912"
CONCEPT_DOI = "10.5281/zenodo.4050912"
PRODUCTION_API = "https://zenodo.org"
SANDBOX_API = "https://sandbox.zenodo.org"

DEFAULT_SAMPLE_ID = "WES-T001-COAD-ART-GL-1"
DEFAULT_ASSEMBLY = "grch38"

REQUIRED_SUFFIXES = [
    "html",
    "xlsx",
    "conf.yaml",
    "log",
    "vcf.gz",
    "vcf.gz.tbi",
    "pass.tsv.gz",
    "classification.tsv.gz",
]
OPTIONAL_SUFFIXES = [
    "biomarker_evidence.tsv.gz",
    "pgx_findings.tsv.gz",
    "secondary_findings.tsv.gz",
]
BINARY_SUFFIXES = ("vcf.gz.tbi", "xlsx", "pdf")
# A home directory path, but not '/home/...' inside a URL such as
# https://www.ncbi.nlm.nih.gov/home/about/
LOCAL_PATH_RE = re.compile(
    r"(?<![\w.\-])(?:/Users/[^/\s]+|/home/[^/\s]+)|[A-Za-z]:\\Users\\"
)
RELEASE_VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")

DESCRIPTION_TEMPLATE = (
    "<ul>\n"
    '<li>Example output files from <a href="https://sigven.github.io/cpsr">'
    "Cancer Predisposition Sequencing Reporter</a></li>\n"
    "<li>Example input variant dataset is artificial - <strong>not</strong> "
    "true germline calls from a patient (hence more findings than what is "
    "usually expected)</li>\n"
    "<li>Generated with CPSR v{version} (genome assembly: {assembly})</li>\n"
    "</ul>"
)


class ReleaseError(Exception):
    """A problem that should stop the release, reported without a traceback."""


# --------------------------------------------------------------------------
# Local checks
# --------------------------------------------------------------------------


def md5sum(path):
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def conf_value(conf_text, pattern):
    match = re.search(pattern, conf_text, flags=re.MULTILINE)
    return match.group(1).strip() if match else None


def xlsx_contains_version(path, version):
    """The xlsx records the R package version that actually built the report."""
    needle = re.compile(r">\s*" + re.escape(version) + r"\s*<")
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if name.endswith(".xml") and needle.search(
                archive.read(name).decode("utf-8", "replace")
            ):
                return True
    return False


def local_path_hits(path):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if LOCAL_PATH_RE.search(line))


def select_files(report_dir, prefix, include_pdf):
    suffixes = REQUIRED_SUFFIXES + OPTIONAL_SUFFIXES
    if include_pdf:
        suffixes = suffixes + ["pdf"]
    return [
        report_dir / f"{prefix}.{suffix}"
        for suffix in suffixes
        if (report_dir / f"{prefix}.{suffix}").is_file()
    ]


def run_checks(args):
    """Validate everything locally and return the files and metadata to send."""
    problems, warnings = [], []
    report_dir = Path(args.report_dir).expanduser().resolve()
    prefix = f"{args.sample_id}.cpsr.{args.assembly}"
    version = args.version

    if not RELEASE_VERSION_RE.match(version):
        problems.append(f"version '{version}' is not a release version (X.Y.Z)")

    description_file = Path(args.package_dir) / "DESCRIPTION"
    package_version = conf_value(
        description_file.read_text(encoding="utf-8"), r"^Version:\s*(\S+)"
    )
    if package_version != version:
        problems.append(
            f"DESCRIPTION has version {package_version}, not {version}; "
            "bump and build the report from the release code first"
        )

    if "-ART-" not in args.sample_id:
        problems.append(
            f"sample '{args.sample_id}' is not marked artificial (-ART-); "
            "only the artificial example may be published"
        )

    for suffix in REQUIRED_SUFFIXES:
        if not (report_dir / f"{prefix}.{suffix}").is_file():
            problems.append(f"missing required file: {prefix}.{suffix}")

    files = select_files(report_dir, prefix, args.include_pdf)
    for path in files:
        if path.stat().st_size == 0:
            problems.append(f"empty file: {path.name}")

    conf_path = report_dir / f"{prefix}.conf.yaml"
    if conf_path.is_file():
        conf_text = conf_path.read_text(encoding="utf-8", errors="replace")
        conf_sample = conf_value(conf_text, r"^sample_id:\s*(.+)$")
        conf_assembly = conf_value(conf_text, r"^genome_assembly:\s*(.+)$")
        conf_version = conf_value(conf_text, r"^\s+cpsr_version:\s*(.+)$")
        if conf_sample != args.sample_id:
            problems.append(f"conf.yaml sample_id is '{conf_sample}'")
        if conf_assembly != args.assembly:
            problems.append(f"conf.yaml genome_assembly is '{conf_assembly}'")
        if conf_version != version:
            warnings.append(
                f"conf.yaml cpsr_version is {conf_version} (set by the Python "
                "wrapper, so it can differ; the xlsx check is authoritative)"
            )

    xlsx_path = report_dir / f"{prefix}.xlsx"
    if xlsx_path.is_file() and not xlsx_contains_version(xlsx_path, version):
        problems.append(f"{xlsx_path.name} does not record CPSR version {version}")

    path_hits = {}
    for path in files:
        if not path.name.endswith(BINARY_SUFFIXES):
            hits = local_path_hits(path)
            if hits:
                path_hits[path.name] = hits
    if path_hits:
        listing = ", ".join(f"{name} ({n})" for name, n in path_hits.items())
        message = f"local paths found in: {listing}"
        if args.allow_local_paths:
            warnings.append(message + " (allowed by --allow-local-paths)")
        else:
            problems.append(
                message + "; generate the report from a neutral path, or pass "
                "--allow-local-paths"
            )

    for warning in warnings:
        print(f"WARNING: {warning}")
    if problems:
        details = "\n".join(f"  - {problem}" for problem in problems)
        raise ReleaseError(f"{len(problems)} check(s) failed:\n{details}")

    metadata = {
        "version": version,
        "publication_date": args.date,
        "description": DESCRIPTION_TEMPLATE.format(
            version=version, assembly=args.assembly.replace("grch", "GRCh")
        ),
    }
    expected = {path.name: md5sum(path) for path in files}
    print(f"Checks passed for {len(files)} files in {report_dir}")
    for path in files:
        print(f"  {path.name}  {path.stat().st_size:,} bytes")
    return files, metadata, expected


# --------------------------------------------------------------------------
# Zenodo API
# --------------------------------------------------------------------------


class SameHostRedirects(urllib.request.HTTPRedirectHandler):
    """Follow redirects only within the API host, so the token stays there."""

    def __init__(self, host):
        self.host = host

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urllib.parse.urlsplit(newurl).netloc != self.host:
            raise ReleaseError(f"refusing redirect to another host: {newurl}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


class Zenodo:
    def __init__(self, base, token):
        self.base = base.rstrip("/")
        self.host = urllib.parse.urlsplit(self.base).netloc
        self.token = token
        self.opener = urllib.request.build_opener(SameHostRedirects(self.host))

    def request(self, method, url, data=None, headers=None, authenticated=True):
        if not url.startswith("http"):
            url = self.base + url
        if urllib.parse.urlsplit(url).netloc != self.host:
            raise ReleaseError(f"refusing to contact another host: {url}")
        request = urllib.request.Request(url, data=data, method=method)
        for key, value in (headers or {}).items():
            request.add_header(key, value)
        if authenticated:
            request.add_header("Authorization", f"Bearer {self.token}")
        try:
            with self.opener.open(request, timeout=600) as response:
                body = response.read()
        except urllib.error.HTTPError as error:
            detail = error.read().decode("utf-8", "replace")[:800]
            raise ReleaseError(f"{method} {url} -> {error.code}: {detail}") from None
        except urllib.error.URLError as error:
            raise ReleaseError(f"{method} {url} failed: {error.reason}") from None
        return json.loads(body) if body else {}

    def send_json(self, method, url, payload):
        return self.request(
            method,
            url,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
        )

    def upload(self, bucket_url, path):
        target = f"{bucket_url}/{urllib.parse.quote(path.name)}"
        with open(path, "rb") as handle:
            return self.request(
                "PUT",
                target,
                data=handle,
                headers={
                    "Content-Type": "application/octet-stream",
                    "Content-Length": str(path.stat().st_size),
                },
            )

    def public_record(self, record_id):
        return self.request(
            "GET", f"/api/records/{record_id}", authenticated=False
        )

    def deposition(self, deposition_id):
        return self.request("GET", f"/api/deposit/depositions/{deposition_id}")


def client(args):
    if args.sandbox:
        base, variable = SANDBOX_API, "ZENODO_SANDBOX_TOKEN"
    else:
        base, variable = PRODUCTION_API, "ZENODO_TOKEN"
    test_base = os.environ.get("ZENODO_TEST_API_BASE")
    if test_base:
        if not test_base.startswith("http://127.0.0.1"):
            raise ReleaseError("ZENODO_TEST_API_BASE must be a local address")
        base = test_base
    token = os.environ.get(variable)
    if not token:
        raise ReleaseError(f"environment variable {variable} is not set")
    print(f"Using {base}")
    return Zenodo(base, token)


def remote_checksums(entries):
    """Map filename to md5 for both deposit and public record file listings."""
    checksums = {}
    for entry in entries:
        name = entry.get("filename") or entry.get("key")
        checksums[name] = entry["checksum"].replace("md5:", "")
    return checksums


def compare(actual_metadata, actual_files, metadata, expected):
    """Return every difference between what Zenodo holds and what was intended."""
    differences = []
    for key, value in metadata.items():
        if actual_metadata.get(key) != value:
            differences.append(
                f"metadata '{key}': expected {value!r}, "
                f"found {actual_metadata.get(key)!r}"
            )
    for name in sorted(set(expected) | set(actual_files)):
        if name not in actual_files:
            differences.append(f"file missing on Zenodo: {name}")
        elif name not in expected:
            differences.append(f"unexpected file on Zenodo: {name}")
        elif expected[name] != actual_files[name]:
            differences.append(f"checksum differs: {name}")
    return differences


def show_metadata(metadata):
    print("Draft metadata:")
    for key in sorted(metadata):
        if key != "prereserve_doi":
            print(f"  {key}: {json.dumps(metadata[key], ensure_ascii=False)}")


def start_draft(api, args):
    """Return (draft id, version of the record it replaces)."""
    if args.new_record:
        if not args.sandbox:
            raise ReleaseError("--new-record is only allowed with --sandbox")
        created = api.send_json(
            "POST",
            "/api/deposit/depositions",
            {
                "metadata": {
                    "title": "Cancer Predisposition Sequencing Reporter (CPSR) "
                    "- example output",
                    "upload_type": "dataset",
                    "creators": [{"name": "Nakken, Sigve"}],
                    "license": "cc-by-4.0",
                    "description": "Sandbox test record",
                }
            },
        )
        return str(created["id"]), None

    record_id = args.record_id or CONCEPT_RECORD_ID
    if args.sandbox and not args.record_id:
        raise ReleaseError("on the sandbox, pass --record-id or --new-record")
    latest = api.public_record(record_id)
    if not args.sandbox and str(latest.get("conceptrecid")) != CONCEPT_RECORD_ID:
        raise ReleaseError(
            f"record {latest.get('id')} is not part of concept {CONCEPT_RECORD_ID}"
        )
    latest_meta = latest.get("metadata", {})
    if latest_meta.get("version") == args.version or (
        f"v{args.version}" in re.sub(r"<[^>]+>", " ", latest_meta.get("description", ""))
    ):
        raise ReleaseError(
            f"latest record {latest['id']} already describes version {args.version}"
        )
    print(f"Latest published record: {latest['id']} ({latest_meta.get('title')})")
    response = api.request(
        "POST", f"/api/deposit/depositions/{latest['id']}/actions/newversion"
    )
    draft_url = response["links"]["latest_draft"]
    return draft_url.rstrip("/").rsplit("/", 1)[-1], latest_meta.get("version")


def cmd_check(args):
    run_checks(args)
    print("No network calls were made.")


def cmd_prepare(args):
    files, metadata, expected = run_checks(args)
    api = client(args)
    draft_id, _ = start_draft(api, args)
    draft = api.deposition(draft_id)
    print(f"Draft {draft_id} created")

    # A new version inherits the previous files, so clear them explicitly.
    for entry in draft.get("files", []):
        api.request(
            "DELETE", f"/api/deposit/depositions/{draft_id}/files/{entry['id']}"
        )
        print(f"  removed inherited {entry['filename']}")

    bucket = draft["links"]["bucket"]
    for path in files:
        result = api.upload(bucket, path)
        remote = result.get("checksum", "").replace("md5:", "")
        if remote != expected[path.name]:
            raise ReleaseError(f"upload checksum mismatch: {path.name}")
        print(f"  uploaded {path.name}")

    # Keep inherited fields such as title, creators and license; replace only
    # the ones that must change with each release.
    merged = dict(draft["metadata"])
    merged.update(metadata)
    api.send_json("PUT", f"/api/deposit/depositions/{draft_id}", {"metadata": merged})

    draft = api.deposition(draft_id)
    differences = compare(
        draft["metadata"], remote_checksums(draft.get("files", [])), metadata, expected
    )
    if differences:
        details = "\n".join(f"  - {item}" for item in differences)
        raise ReleaseError(f"draft {draft_id} does not match:\n{details}")

    show_metadata(draft["metadata"])
    print("\nDraft verified. Nothing is public yet.")
    print(f"Review: {draft['links'].get('html', '(no link returned)')}")
    print(
        f"Then:   python3 data-raw/zenodo_release.py publish --draft-id {draft_id} "
        f"--version {args.version} --report-dir {args.report_dir}"
        + (" --sandbox" if args.sandbox else "")
    )


def cmd_publish(args):
    _, metadata, expected = run_checks(args)
    api = client(args)
    draft = api.deposition(args.draft_id)
    if draft.get("submitted") and draft.get("state") == "done":
        raise ReleaseError(f"deposition {args.draft_id} is already published")

    # The draft's date is kept, so a review on a later day doesn't alter it.
    metadata["publication_date"] = draft["metadata"].get(
        "publication_date", metadata["publication_date"]
    )
    differences = compare(
        draft["metadata"], remote_checksums(draft.get("files", [])), metadata, expected
    )
    if differences:
        details = "\n".join(f"  - {item}" for item in differences)
        raise ReleaseError(f"draft does not match local files:\n{details}")

    show_metadata(draft["metadata"])
    target = "SANDBOX" if args.sandbox else "PRODUCTION"
    print(f"\nThis publishes {len(expected)} files to {target} and mints a DOI.")
    print("Published files can't be changed or replaced.")
    phrase = f"publish {args.version}"
    try:
        answer = input(f"Type '{phrase}' to continue: ")
    except EOFError:
        answer = ""
    if answer.strip() != phrase:
        raise ReleaseError("not confirmed; nothing was published")

    published = api.request(
        "POST", f"/api/deposit/depositions/{args.draft_id}/actions/publish"
    )
    record_id = published.get("record_id") or published["id"]
    print(f"Published record {record_id}, DOI {published.get('doi')}")

    # Public indexing can lag briefly, so retry before reporting a difference.
    differences = ["record not checked yet"]
    for attempt in range(args.verify_attempts):
        if attempt:
            time.sleep(args.verify_wait)
        try:
            record = api.public_record(record_id)
        except ReleaseError as error:
            differences = [str(error)]
            continue
        differences = compare(
            record.get("metadata", {}),
            remote_checksums(record.get("files", [])),
            metadata,
            expected,
        )
        if not differences:
            break
    if differences:
        details = "\n".join(f"  - {item}" for item in differences)
        raise ReleaseError(
            f"published, but the public record doesn't match yet:\n{details}\n"
            f"Check {api.base}/records/{record_id}"
        )
    print("Public record verified: metadata and all checksums match.")

    if not args.sandbox:
        latest = api.public_record(CONCEPT_RECORD_ID)
        if str(latest.get("id")) == str(record_id):
            print(f"Concept DOI {CONCEPT_DOI} resolves to this version.")
        else:
            print(
                f"WARNING: concept record still resolves to {latest.get('id')}; "
                "this can lag, so check again shortly."
            )


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    today = datetime.date.today().isoformat()
    package_dir = Path(__file__).resolve().parent.parent

    for name, handler in (
        ("check", cmd_check),
        ("prepare", cmd_prepare),
        ("publish", cmd_publish),
    ):
        sub = commands.add_parser(name)
        sub.set_defaults(handler=handler)
        sub.add_argument("--version", required=True, help="release version, X.Y.Z")
        sub.add_argument("--report-dir", required=True)
        sub.add_argument("--sample-id", default=DEFAULT_SAMPLE_ID)
        sub.add_argument("--assembly", default=DEFAULT_ASSEMBLY)
        sub.add_argument("--date", default=today, help="publication date, YYYY-MM-DD")
        sub.add_argument("--package-dir", default=str(package_dir))
        sub.add_argument("--include-pdf", action="store_true")
        sub.add_argument(
            "--allow-local-paths",
            action="store_true",
            help="publish even if files contain paths such as /Users/<name>",
        )
        if name != "check":
            sub.add_argument("--sandbox", action="store_true")
        if name == "prepare":
            sub.add_argument(
                "--record-id", help="record to version; defaults to the CPSR concept"
            )
            sub.add_argument(
                "--new-record",
                action="store_true",
                help="sandbox only: create a fresh test record",
            )
        if name == "publish":
            sub.add_argument("--draft-id", required=True)
            sub.add_argument("--verify-attempts", type=int, default=6)
            sub.add_argument("--verify-wait", type=float, default=10.0)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        args.handler(args)
    except ReleaseError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
