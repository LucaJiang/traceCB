# Contributing to traceCB

Please open an issue before making a substantial change to the statistical
method, input schema, or paper analysis workflow. Keep scientific behavior
changes separate from formatting-only changes so that reviewers can audit the
effect on results.

## Development setup

```bash
python -m pip install -e ".[ci,docs]"
```

Before submitting a pull request, run:

```bash
python -m pytest -q
python -m build
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
