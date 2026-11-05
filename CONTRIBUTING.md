# Contributing

Thank you for looking at PeatGuard NISAR. Contributions that improve the
radar processing, the scientific audit trail or the viewer are welcome.

## Ground rules

- Every method or physical claim in code comments, docs or commit messages
  cites the named paper it comes from and carries a confidence tag:
  `[E]` established, `[S]` standard inference, `[H]` hypothesis. Do not
  conflate sources (Tay et al. 2025 is radar subsidence; Dadap et al. 2021 is
  an optical canal map).
- Never present an InSAR product as a water-table proxy, and never calibrate
  our estimates to the Tay et al. 2025 benchmark; it is an independent check.
- Keep the maturity labels (`PROVISIONAL`, `EXPERIMENTAL`, `not_ground_truth`)
  on every NISAR raster and catalogue entry.
- No credentials in the tree. Earthdata, CDS and other tokens are read from
  `~/.netrc`, `~/.cdsapirc` or environment variables.
- Plain ASCII punctuation in prose: no em dashes, no emojis.

## Workflow

1. Fork and branch from `main`.
2. `make install` (or `pip install -e '.[ci,dashboard]'`), then `make check`.
   The suite runs offline; it needs no satellite data or cluster.
3. Add or update tests under `tests/` (or `cluster/full/tests/` for cluster
   job regressions) for any behaviour change.
4. Open a pull request with a short description of what changed, why, and
   which document or paper backs any scientific change.

Commits use plain `git commit -m` messages; no generated co-author trailers.

## Reporting problems

Open an issue with the command you ran, the Python version, the traceback,
and, for data issues, the GUNW granule name. Do not attach raw NISAR products.
