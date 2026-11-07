"""Local artifact preparation; publication is deliberately handled separately."""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import tomllib
import zipfile
from pathlib import Path

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.version import Version

ROOT = Path(__file__).resolve().parents[1]


def read_toml(path):
    with path.open("rb") as stream:
        return tomllib.load(stream)


def run(*args, cwd=ROOT):
    subprocess.run([str(arg) for arg in args], cwd=cwd, check=True)


def conda_python_dependencies(project, mapping):
    overrides = mapping.get("overrides", {})
    reasons = mapping.get("override_reasons", {})
    unknown = set(overrides) - set(project["dependencies"])
    if unknown:
        raise ValueError(f"Conda overrides have no Python requirement: {sorted(unknown)}")
    if set(overrides) != set(reasons):
        raise ValueError("Every Conda override must have exactly one documented reason")
    if any(not reason.strip() for reason in reasons.values()):
        raise ValueError("Conda override reasons must not be empty")

    translated = []
    for raw_requirement in project["dependencies"]:
        requirement = Requirement(raw_requirement)
        if raw_requirement in overrides:
            replacements = overrides[raw_requirement]
            if not replacements:
                raise ValueError(f"Empty Conda override for {raw_requirement}")
            first = Requirement(replacements[0].replace(" ", ""))
            if first.name.lower().replace("_", "-") != requirement.name.lower().replace("_", "-"):
                raise ValueError(f"Conda override must retain package {requirement.name}")
            for replacement in replacements:
                Requirement(replacement.replace(" ", ""))
            translated.extend(replacements)
            continue

        if requirement.extras or requirement.marker or requirement.url:
            raise ValueError(
                f"Python requirement needs an explicit Conda override: {raw_requirement}"
            )
        specifiers = list(requirement.specifier)
        if any(spec.operator != ">=" for spec in specifiers):
            raise ValueError(
                f"Python version syntax needs an explicit Conda override: {raw_requirement}"
            )
        suffix = f" {requirement.specifier}" if requirement.specifier else ""
        translated.append(f"{requirement.name}{suffix}")
    return translated


def check(root=ROOT):
    project = read_toml(root / "pyproject.toml")["project"]
    if Version("3.11") not in SpecifierSet(project["requires-python"]):
        raise ValueError("Update the Conda host/test Python target for this project's Python range")
    mapping = read_toml(ROOT / "conda-recipe/dependencies.toml")
    conda_python_dependencies(project, mapping)
    if project["scripts"] != {"gmri2fem": "_cli:gmri2fem"}:
        raise ValueError("CLI changed: update the Conda entry point and artifact tests")
    policy = read_toml(ROOT / "conda-recipe/fsl.toml")
    if Version(project["version"]) not in SpecifierSet(policy["package"]["gmri2fem"]):
        raise ValueError("Companion compatibility range excludes this gmri2fem version")
    return project


def stage(sdist, directory):
    if directory.exists():
        raise FileExistsError(f"Refusing to overwrite {directory}; choose a fresh staging directory")
    directory.mkdir(parents=True)
    with tarfile.open(sdist) as archive:
        members = archive.getmembers()
        roots = {Path(member.name).parts[0] for member in members}
        if len(roots) != 1:
            raise ValueError("Expected one top-level sdist directory")
        # The data filter rejects traversal paths and unsafe links.
        archive.extractall(directory, filter="data")
    (directory / roots.pop()).rename(directory / "source")
    check(directory / "source")
    shutil.copytree(ROOT / "conda-recipe", directory / "conda-recipe")
    return directory / "source"


def checksums(directory):
    artifacts = sorted(
        path for path in directory.rglob("*")
        if path.is_file() and (path.suffix in {".whl", ".conda"} or path.name.endswith((".tar.gz", ".tar.bz2")))
    )
    (directory / "SHA256SUMS").write_text("".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(directory)}\n"
        for path in artifacts
    ))


def build():
    project = check()
    dist = ROOT / "dist"
    if dist.exists() and any(dist.iterdir()):
        raise FileExistsError("dist is not empty; preserve old artifacts and use a clean checkout")
    run(sys.executable, "-m", "build", "--sdist")
    sdist, = dist.glob("*.tar.gz")
    source = stage(sdist, ROOT / "build/release")
    run(sys.executable, "-m", "build", "--wheel", "--outdir", dist, cwd=source)
    wheel, = dist.glob("*.whl")
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        required = {"_cli.py"} | {f"{name}/__init__.py" for name in ("gmri2fem", "i2m", "dti", "brainmeshing")}
        if not required <= names:
            raise ValueError(f"Missing wheel payload: {required - names}")
    run(sys.executable, "-m", "twine", "check", *dist.glob("*.whl"), sdist)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=ROOT, text=True,
    ).strip())
    (dist / "provenance.json").write_text(json.dumps({
        "version": project["version"], "commit": commit,
        "dirty_worktree": dirty,
        "sdist_sha256": hashlib.sha256(sdist.read_bytes()).hexdigest(),
    }, indent=2) + "\n")
    checksums(dist)


def test_conda(work, output, companion_only):
    project = read_toml(work / "source/pyproject.toml")["project"]
    policy = read_toml(work / "conda-recipe/fsl.toml")
    for kind in (["fsl"] if companion_only else ["native", "fsl"]):
        directory = work / f"test-{kind}"
        directory.mkdir(exist_ok=True)
        channels = [
            output.as_uri(), "conda-forge/label/fenics-dev", "conda-forge",
            "jorgenriseth", "bioconda",
        ]
        if kind == "fsl":
            channels.append("https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public")
        name = policy["package"]["name"] if kind == "fsl" else project["name"]
        version = policy["package"]["version"] if kind == "fsl" else project["version"]
        manifest = directory / "pixi.toml"
        manifest.write_text(
            '[workspace]\nname = "artifact-test"\nplatforms = ["linux-64"]\n'
            f"channels = {json.dumps(channels)}\n"
            '[dependencies]\npython = "3.11.*"\npip = "*"\npytest = ">=8,<9"\n'
            f'{name} = {{ version = "=={version}", channel = "{output.as_uri()}" }}\n'
            + (
                '[activation.env]\nFSLDIR = "$CONDA_PREFIX"\nFSLOUTPUTTYPE = "NIFTI_GZ"\n'
                if kind == "fsl" else ""
            )
        )
        run("pixi", "install", "--manifest-path", manifest, "--run-post-link-scripts")
        prefix = ("pixi", "run", "--manifest-path", manifest, "--locked")
        run(*prefix, "python", "-m", "pip", "check", cwd=directory)
        run(*prefix, "python", "-c",
            "from importlib.metadata import version; print('Installed gmri2fem:', version('gmri2fem'))",
            cwd=directory)
        for command in [
            ["--help"], ["i2m", "concentrations2mesh", "--help"],
            ["i2m", "dti2mesh", "--help"], ["i2m", "collect", "--help"],
            ["brainmeshing", "--help"],
        ]:
            run(*prefix, "gmri2fem", *command, cwd=directory)
        tests = ["test_native_integration_smoke.py"]
        if kind == "fsl":
            tests.append("test_fsl.py")
        run(*prefix, "python", "-m", "pytest", "-q", "--import-mode=importlib",
            *(work / "source/tests" / test for test in tests), cwd=directory)
        record = ROOT / "dist/test-environments" / kind
        record.mkdir(parents=True, exist_ok=True)
        for filename in ("pixi.toml", "pixi.lock"):
            shutil.copyfile(directory / filename, record / filename)


def conda(sdist, companion_only=False, directory=None):
    work = directory or ROOT / "build/conda-release"
    stage(sdist, work)
    output = ROOT / "dist/conda"
    channels = [
        "-c", output, "-c", "conda-forge/label/fenics-dev", "-c", "conda-forge",
        "-c", "jorgenriseth", "-c", "bioconda",
    ]
    output.mkdir(parents=True, exist_ok=True)
    if not companion_only:
        run("conda", "build", work / "conda-recipe", "--no-anaconda-upload",
            "--no-test", "--python", "3.11", "--override-channels", *channels[2:],
            "--output-folder", output)
    run("conda", "index", output)
    run("conda", "build", work / "conda-recipe/fsl", "--no-anaconda-upload",
        "--no-test", "--override-channels", *channels,
        "-c", "https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public",
        "--output-folder", output)
    run("conda", "index", output)
    test_conda(work, output, companion_only)
    checksums(ROOT / "dist")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["check", "build", "conda", "stage", "test-conda"])
    parser.add_argument("--sdist", type=Path)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--companion-only", action="store_true")
    parser.add_argument("--tag", help="Require a library release tag matching project.version")
    args = parser.parse_args()
    project = check()
    if args.tag and args.tag != f"v{project['version']}":
        parser.error(f"Expected tag v{project['version']}, got {args.tag}")
    if args.command == "build":
        build()
    elif args.command == "test-conda":
        if args.directory is None:
            parser.error("--directory is required")
        test_conda(args.directory.resolve(), ROOT / "dist/conda", args.companion_only)
        checksums(ROOT / "dist")
    elif args.command in {"conda", "stage"}:
        if args.sdist is None:
            parser.error("--sdist is required")
        if args.command == "stage":
            if args.directory is None:
                parser.error("--directory is required")
            stage(args.sdist.resolve(), args.directory.resolve())
        else:
            conda(args.sdist.resolve(), args.companion_only,
                  args.directory.resolve() if args.directory else None)


if __name__ == "__main__":
    main()
