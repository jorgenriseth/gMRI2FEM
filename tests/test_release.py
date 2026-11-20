import importlib.util
from pathlib import Path
import tarfile

from jinja2 import Environment, StrictUndefined
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("release", ROOT / "scripts/release.py")
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


def test_metadata_consistency():
    assert release.check()["name"] == "gmri2fem"


@pytest.mark.parametrize("recipe", ["meta.yaml", "fsl/meta.yaml"])
def test_recipe_structure(recipe):
    recipe_path = ROOT / "conda-recipe" / recipe

    def load_file_data(path):
        path = Path(path)
        if path.name == "pyproject.toml":
            path = ROOT / "pyproject.toml"
        return release.read_toml(path)

    template = Environment(undefined=StrictUndefined).from_string(recipe_path.read_text())
    rendered = template.render(
        RECIPE_DIR=str(recipe_path.parent), load_file_data=load_file_data, PYTHON="python",
    )
    metadata = yaml.safe_load(rendered)
    assert "test" in metadata
    assert "test" not in metadata["requirements"]
    assert metadata["test"]["commands"]
    assert "license_file" in metadata["about"]
    if recipe == "meta.yaml":
        project = release.check()
        mapping = release.read_toml(ROOT / "conda-recipe/dependencies.toml")
        assert metadata["package"]["name"] == project["name"]
        assert metadata["package"]["version"] == project["version"]
        assert metadata["build"]["entry_points"] == ["gmri2fem = _cli:gmri2fem"]
        assert f"python {project['requires-python']}" in metadata["requirements"]["run"]
        run_requirements = metadata["requirements"]["run"][1:]
        expected = release.conda_python_dependencies(project, mapping)
        expected.extend(
            f"{name} {constraint}" for name, constraint in mapping["native"].items()
        )
        assert run_requirements == expected


def test_stage_rejects_existing_directory(tmp_path):
    with pytest.raises(FileExistsError):
        release.stage(tmp_path / "unused.tar.gz", tmp_path)


def test_stage_rejects_traversal(tmp_path):
    path = tmp_path / "bad.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo("../outside")
        archive.addfile(member)
    with pytest.raises(tarfile.FilterError):
        release.stage(path, tmp_path / "stage")
