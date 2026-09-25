# Tests

Run the tracked automated tests from the repository root:

```bash
python -m pytest -q
```

`test_core.py` covers package metadata and numerical safety checks for the core
estimators. `test_esnp_direction_replicate.py` exercises allele alignment and
sign handling in `src/figures/esnp_replication.py`. The tutorial notebook is also
executed in CI as an end-to-end smoke test. The remaining tests cover
chromosome 22 CLI/plotting behavior, simulation variance terms, toy-data LDSC,
and independence of repeated LDSC fits. CI runs the tests from the source
archive against the installed wheel.

Long-running simulation grids and their generated outputs are intentionally not
part of the unit-test suite; use the documented entry points under
`src/simulation/` for those analyses.
