# Tests

Run the tracked automated tests from the repository root:

```bash
python -m pytest -q
```

`test_core.py` covers package metadata and numerical safety checks for the core
estimators. `test_esnp_direction_replicate.py` exercises allele alignment and
sign handling in `src/figures/esnp_replication.py`. The tutorial notebook is also
executed in CI as an end-to-end smoke test.

Long-running simulation grids and their generated outputs are intentionally not
part of the unit-test suite; use the documented entry points under
`src/simulation/` for those analyses.
