# Contributing to CPSR

Thanks for your interest in contributing to the Cancer Predisposition Sequencing
Reporter (CPSR). Contributions of all kinds are welcome — bug reports,
documentation improvements, feature suggestions, and code.

This document explains how to report problems and how to propose changes.

## Code of conduct

This project follows a [Code of Conduct](CODE_OF_CONDUCT.md). By participating,
you are expected to uphold it. Please report unacceptable behaviour to the
maintainers.

## Ways to contribute

- **Report a bug** — open an issue (see below for what to include)
- **Request a feature** — open an issue describing the use case
- **Ask a question** — use GitHub Discussions or open an issue
- **Improve documentation** — corrections and clarifications are very welcome
- **Contribute code** — see *Development workflow* below

## Reporting bugs

Most issues we receive can only be diagnosed with the full run context, so
please include **all** of the following. Issues without this information will
usually need a follow-up before we can help.

1. **CPSR version**
2. **Reference data bundle version** and the genome assembly used
   (`grch37` or `grch38`)
3. **Installation method** — conda/mamba, Docker, Singularity/Apptainer —
   including the image tag or environment specification
4. **Operating system** and, if relevant, the compute environment
   (laptop, HPC/cluster, cloud)
5. **The exact command you ran**, in full, with all arguments
6. **The complete error message and log output** (please paste as text in a
   fenced code block rather than as a screenshot)
7. **Input and configuration details** — variant caller used, the virtual gene
   panel selected (or a custom panel/BED, if used), whether secondary findings
   and/or GWAS variants were enabled, and the approximate number of variants
8. **A minimal example input** that reproduces the problem, if you are able to
   share one

Please do not attach patient-identifiable data to a public issue. Germline
variant data is particularly sensitive; if a problem can only be reproduced with
such data, say so in the issue and we will find another way to investigate.

## Requesting features

When proposing a feature, describe the clinical or research use case rather than
only the implementation. CPSR classifies germline variants against ACMG/AMP
criteria using public resources (ClinVar, gnomAD, cancerhotspots.org), so proposals that
affect variant classification or evidence weighting should
reference the relevant guideline, gene–disease validity source, or evidence
resource.

## Development workflow

We use a simple two-branch model:

- **`main`** — reflects the current released version. Do not open pull requests
  against `main`.
- **`dev`** — the active development branch. **All pull requests should target
  `dev`.**

To contribute code:

1. Fork the repository and create a branch from `dev`
   (e.g. `feature/short-description` or `fix/short-description`)
2. Make your changes, keeping the pull request focused on a single concern
3. Update documentation and the changelog where relevant
4. Open a pull request against `dev`, describing **what** the change does and
   **why** it is needed

For anything substantial — new evidence sources, changes to ACMG/AMP
classification logic, changes to panel handling, or architectural changes —
please open an issue to discuss the approach before writing code. This saves
time on both sides.

## Code organisation

CPSR is built on the same framework as
[PCGR](https://github.com/sigven/pcgr) and shares much of its infrastructure:

- **Python** — the pipeline, variant annotation, and orchestration of the run
- **R** (the `cpsr` package) — report generation and visualisation via Quarto,
  together with the ACMG/AMP-based variant classification logic. Some shared
  functionality is inherited from the [`pcgrr`](https://github.com/sigven/pcgr) R
  package, which CPSR builds on.

Because classification and reporting both live in the `cpsr` R package, changes
to how variants are classified and changes to how results are presented can
touch the same code. If you are unsure where a change belongs, open an issue and
ask before starting work.

Note that CPSR is distributed together with PCGR, so some changes may affect
both tools. Please say so in your pull request if you believe this is the case.

## Coding conventions

- Follow the style of the surrounding code
- R: follow the existing style in `cpsr`
- Python: follow PEP 8 where practical; add docstrings to new functions
- Keep commits reasonably self-contained with descriptive messages
- Do **not** commit reference data bundles, large binary files, or test data
  containing patient information

## Testing and validation

Automated test coverage is currently limited and is actively being expanded. In
the meantime, we ask contributors to describe how a change was validated:

- State the command(s) you ran to exercise the change
- For changes affecting annotation or classification, show the effect on output
  — for example, the relevant variants or classification fields before and after
- Confirm that an end-to-end run completes on example data
- Where tests exist, please run them and add new ones covering your change

Continuous integration runs automatically on pull requests. Please make sure it
passes before requesting review.

## Releases and versioning

CPSR follows semantic versioning. Releases are cut from `dev` into `main` by the
maintainers, accompanied by an updated changelog. Reference data bundles are
versioned separately and tied to specific releases — see the documentation for
the compatible bundle for each version.

## Licence

CPSR is released under the MIT licence. By contributing, you agree that your
contributions will be licensed under the same terms.

## Getting in touch

- **Bugs and feature requests** — GitHub Issues
- **Questions and general discussion** — GitHub Discussions (to come)
- **Documentation** — https://sigven.github.io/cpsr/

Thanks for your effort to improve CPSR!
