# Contributing to traceCB

Please open an issue before making a substantial change to the statistical
method, input schema, or paper analysis workflow. Keep scientific behavior
changes separate from formatting-only changes so that reviewers can audit the
effect on results.

## Development setup

```bash
python -m pip install -e ".[ci,test,docs]"
```

Before submitting a pull request, run:

```bash
python -m pytest -q
python -m build
python -m twine check --strict dist/*
mkdocs build --strict
jupyter nbconvert --to notebook --execute \
  docs/tutorial/run_traceCB.ipynb \
  --output run_traceCB_checked.ipynb \
  --output-dir /tmp \
  --ExecutePreprocessor.timeout=600
```

Do not commit controlled-access data, credentials, local absolute-path
configuration, generated result directories, or notebook execution caches.
When an analysis result changes, record the command, parameters, random seed,
software revision, and affected manuscript figure or table in the pull request.


## Release artifacts

Build from a clean checkout of the intended release commit. The wheel contains
only the installable `traceCB` library and its metadata. The source archive also
contains workflows, tests, documentation, and the public toy inputs so that it
can be tested independently of Git. Run the tests from an extracted source
archive against the installed wheel, and execute both tutorial notebooks.
CI checks this on Python 3.10 and 3.12; notebook execution uses Python 3.12.

Before publishing, ensure that the release tag points to the reviewed commit
and that `pyproject.toml`, `traceCB.__version__`, and `CITATION.cff` agree. Record
the commit and environment used for manuscript results. Passing the toy example
and unit tests does not replace reproducing the full analyses with their
external datasets and R/PLINK/LDSC dependencies.
