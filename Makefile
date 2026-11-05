# PeatGuard NISAR: developer entry points.
#
#   make install   editable install with the dashboard and test extras
#   make demo      start the interactive COG viewer on data/sample (port 8080)
#   make test      run the offline test suite (no satellite data required)
#   make lint      ruff + syntax check of Python and cluster shell jobs
#
# PY selects the interpreter (default: python3). PYTHONPATH=src lets every
# target run from a clean checkout without installing the package.

PY      ?= python3
PORT    ?= 8080
HOST    ?= 0.0.0.0
SAMPLE  := $(CURDIR)/data/sample
export PYTHONPATH := $(CURDIR)/src:$(PYTHONPATH)
export MPLBACKEND := Agg

.PHONY: install demo test lint check clean

install:
	$(PY) -m pip install -e '.[ci,dashboard]'

demo:
	@echo "PeatGuard NISAR viewer on http://localhost:$(PORT)/  (sample data: $(SAMPLE))"
	PEATGUARD_WEBDATA_DIR=$(SAMPLE) $(PY) -m uvicorn peatguard.dashboard.app:app --host $(HOST) --port $(PORT)

test:
	$(PY) -m pytest -q

lint:
	$(PY) -m ruff check src tests scripts viz || true
	$(PY) scripts/ci/check_pipeline.py

check: lint test

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	rm -rf .pytest_cache build dist src/*.egg-info
