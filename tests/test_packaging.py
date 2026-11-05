"""Keep the offline test environment and installed dashboard resources explicit."""

from __future__ import annotations

import fnmatch
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 remains supported by the package.
    import tomli as tomllib


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _metadata() -> dict:
    with (PROJECT_ROOT / "pyproject.toml").open("rb") as stream:
        return tomllib.load(stream)


def test_ci_extra_declares_test_client_and_avoids_cluster_processing_install() -> None:
    extras = _metadata()["project"]["optional-dependencies"]
    ci_dependencies = extras["ci"]
    assert any(requirement.startswith("pytest>=") for requirement in ci_dependencies)
    assert any(requirement.startswith("httpx>=") for requirement in ci_dependencies)
    assert any(requirement.startswith("build>=") for requirement in ci_dependencies)
    assert any(requirement.startswith("tomli>=") for requirement in ci_dependencies)
    assert all("mintpy" not in requirement.lower() for requirement in ci_dependencies)
    assert all("isce" not in requirement.lower() for requirement in ci_dependencies)
    assert any(requirement.startswith("fastapi>=") for requirement in extras["dashboard"])
    assert any(requirement.startswith("titiler.core>=") for requirement in extras["dashboard"])


def test_required_installed_package_resources_are_declared() -> None:
    patterns = _metadata()["tool"]["setuptools"]["package-data"]["peatguard"]
    for relative in ("data/default.yaml", "dashboard/templates/index.html"):
        assert (PROJECT_ROOT / "src/peatguard" / relative).is_file()
        assert any(fnmatch.fnmatch(relative, pattern) for pattern in patterns), relative


def test_default_pytest_collection_includes_cluster_regressions() -> None:
    paths = _metadata()["tool"]["pytest"]["ini_options"]["testpaths"]
    assert "tests" in paths
    assert "cluster/full/tests" in paths
