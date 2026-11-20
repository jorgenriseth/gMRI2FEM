import importlib

import click
import pytest
from click.testing import CliRunner

from _cli import LazyGroup, gmri2fem
from gmri2fem.external import require_commands


@pytest.mark.parametrize("args", [[], ["mri"], ["seg"], ["dti"], ["i2m"], ["brainmeshing"]])
def test_help_does_not_import_commands(monkeypatch, args):
    original = importlib.import_module

    def guarded(name, *extra, **kwargs):
        assert name not in {"brainmeshing.mesh_generation", "i2m.evaluation_data", "dolfin", "SVMTK"}
        return original(name, *extra, **kwargs)

    monkeypatch.setattr(importlib, "import_module", guarded)
    result = CliRunner().invoke(gmri2fem, [*args, "--help"])
    assert result.exit_code == 0, result.output


def test_lazy_help_never_loads_implementation(monkeypatch):
    group = LazyGroup("example", lazy_subcommands={"native": "nonexistent.command"})
    result = CliRunner().invoke(group, ["--help"])
    assert result.exit_code == 0
    assert "native" in result.output


def test_missing_external_commands(monkeypatch):
    monkeypatch.setattr("gmri2fem.external.shutil.which", lambda _: None)
    with pytest.raises(click.ClickException, match="fslroi.*greedy"):
        require_commands("fslroi", "greedy")


def test_existing_external_commands(monkeypatch):
    monkeypatch.setattr("gmri2fem.external.shutil.which", lambda name: f"/bin/{name}")
    require_commands("fslroi")


def test_missing_native_dependency_is_actionable(monkeypatch):
    def missing(name):
        raise ModuleNotFoundError("No module named dolfin", name="dolfin")

    monkeypatch.setattr(importlib, "import_module", missing)
    result = CliRunner().invoke(gmri2fem, ["i2m", "collect"])
    assert result.exit_code == 1
    assert "Conda/Pixi" in result.output


def test_unrelated_import_errors_are_not_hidden(monkeypatch):
    def missing(name):
        raise ModuleNotFoundError("unexpected import", name="unrelated")

    monkeypatch.setattr(importlib, "import_module", missing)
    result = CliRunner().invoke(gmri2fem, ["i2m", "collect"])
    assert isinstance(result.exception, ModuleNotFoundError)
