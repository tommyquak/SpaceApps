# Processing environment

| component | version | note |
|---|---|---|
| ISCE2 | **2.6.5** | Tay et al. 2025 used 2.6.3 |
| MintPy | **1.6.4** | Tay et al. 2025 used 1.5.3 |
| snaphu | **2.0.7** | binary, NOT the python wrapper alone |
| GDAL | 3.10.3 | |
| numpy | 1.26.4 | |
| h5py | 3.16.0 | |
| python | 3.11.16 | |

`insar-env.lock` is the explicit package list; `bootstrap_insar_env.sh` recreates
the environment and performs the one manual step the lockfile cannot express.

## The step that is not in the lockfile

`conda install snaphu` installs a **python wrapper** and puts **no `snaphu`
command on PATH**. ISCE gates on `shutil.which('snaphu')`; when it is absent,
`Alos2Proc/runUnwrapSnaphu.py` silently calls
`snaphuUnwrap(..., costMode='SMOOTH', initMethod='MCF', initOnly=True)`, running
minimum-cost-flow **initialisation only** with the statistical cost solver never
invoked. That is not an error, it is a fallback, and it produced every
interferogram in this project before 2026-09-03.

The real Stanford binary ships inside the python package at
`site-packages/snaphu/snaphu` and must be symlinked into `bin/`. Verify with
`shutil.which`, not `which`, because `shutil.which` is the test ISCE performs.

## Run it on a compute node

The login node has `ulimit -v` of 1 GB and the solver dies with
`Out of memory allocating ...` partway through. Use `srun`.
