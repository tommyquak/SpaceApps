# Legacy deployment files

`cloud/` holds the original Google Cloud deployment of the pipeline (Cloud
Build image, Cloud Run jobs, a Workflows definition and a Cloud Scheduler
setup script). PeatGuard no longer runs on that project: the pipeline moved to
a Slurm cluster (`cluster/`) and the web stack to a plain Docker Compose host.
The files are kept for reference only.

- They reference a GCP project, region and bucket that are not maintained.
- The scheduler script expects to be run from the old monorepo layout
  (`peatguard/peatguard/`), which is now the repository root.
- The dashboard container definition that is still useful was moved to the
  repository root as `Dockerfile.dashboard` and is wired into
  `docker-compose.yaml`.
- `peatguard.export.gcs` and the optional `cloud` extra in `pyproject.toml`
  remain functional for anyone who wants to sync COGs to a bucket, but nothing
  in the default workflow needs them.
