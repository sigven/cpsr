# Zenodo release steps

Publishes the example report as a new version of the CPSR record on Zenodo,
using `data-raw/zenodo_release.py`. The README and site badges use the concept
DOI (`10.5281/zenodo.4050912`), so they don't need changing for a release.

## Tokens

Set these in your shell. The two sites have separate accounts, and both tokens
need the `deposit:write` and `deposit:actions` scopes.

| Use         | Variable               | Create it at         |
|-------------|------------------------|----------------------|
| Testing     | `ZENODO_SANDBOX_TOKEN` | `sandbox.zenodo.org` |
| Real upload | `ZENODO_TOKEN`         | `zenodo.org`         |

## Steps

Replace `2.4.0` with the version being released.

1. Bump to the release version. The script checks `DESCRIPTION`.

2. Build the example report (`WES-T001-COAD-ART-GL-1`) with the release code.
   The script refuses a report whose `.xlsx` doesn't record this version. A
   Docker run should avoid `/Users/...` paths in the output (not yet tried).

3. Run the local checks:

   ```bash
   python3 data-raw/zenodo_release.py check --version 2.4.0 --report-dir <dir>
   ```

4. Test on the sandbox. The script hasn't run against the real API yet, so do
   this at least once before the first real release. `--new-record` is for a
   sandbox with no copy of the record.

   ```bash
   python3 data-raw/zenodo_release.py prepare --sandbox --new-record --version 2.4.0 --report-dir <dir>
   python3 data-raw/zenodo_release.py publish --sandbox --draft-id <id> --version 2.4.0 --report-dir <dir>
   ```

5. Create the real draft. This prints a review link and the draft ID, and
   nothing is public yet.

   ```bash
   python3 data-raw/zenodo_release.py prepare --version 2.4.0 --report-dir <dir>
   ```

6. Review the draft, then publish. You'll be asked to type `publish 2.4.0`.
   Published files can't be replaced, only superseded by a new version.

   ```bash
   python3 data-raw/zenodo_release.py publish --draft-id <id> --version 2.4.0 --report-dir <dir>
   ```

## Options

- `--allow-local-paths`: publish even if files contain paths like `/Users/<name>`.
- `--include-pdf`: also upload the PDF report, which is left out by default.
