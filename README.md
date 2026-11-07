# gMRI2FEM

`gmri2fem` is a Python package and command-line interface for processing
contrast-enhanced human brain MRI, estimating tracer concentrations, and
mapping image data to finite-element meshes for mathematical modelling of
brain tracer transport. The package was developed for MRI processing in the
GRIP project.

## Installation

Choose the installation route based on the features you need.

### Python-only functionality

For concentration calculations and other Python functionality that does not
need the native meshing stack or FSL:

```bash
python -m pip install gmri2fem
```

This installs the Python package and CLI, but not legacy FEniCS, SVMTK, FSL,
Greedy, or FreeSurfer. Python-only releases target Python 3.11 and 3.12.

### Native meshing and FSL features

Conda recipes for gmri2fem and its optional FSL companion are being prepared
for publication. Until a validated release is available, this installation
route is not yet available from the configured channels.

Once published, use Pixi/Conda for features requiring legacy FEniCS
(`dolfin`, not FEniCSx/`dolfinx`) and SVMTK. Install `gmri2fem-fsl` as well
when you need FSL commands; it will install a compatible gmri2fem version.
Greedy must be installed separately. FreeSurfer/FastSurfer setup and licenses
are also managed separately.

Managed native environments currently target Linux x86-64. A package version
range does not guarantee numerical equivalence across different native stacks.

### Gonzo pipeline

Gonzo is a separate downstream MRI processing pipeline with its own
environment and external-tool setup. See the
[Gonzo repository](https://github.com/jorgenriseth/gonzo) and its documentation
for installation and pipeline instructions; do not use this package's
Python-only pip installation as a substitute for the Gonzo environment.

## Usage

The CLI has nested subcommands. Start by listing available commands:

```bash
gmri2fem --help
```

Example output:

```text
Usage: gmri2fem [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  brainmeshing
  dcm2nii-ll
  dcm2nii-mixed
  dti
  i2m
  mri
  seg
```

To explore MRI commands and their options:

```bash
gmri2fem mri --help
gmri2fem mri looklocker-t1map --help
```

For example, `gmri2fem mri looklocker-t1map --help` describes the input,
timestamps, and output options for that command.

## Contributors and release maintainers

Use Pixi **0.83.0** for development and release tasks:

```bash
pixi install --locked
pixi run --locked test
```

For release tooling and a local artifact rehearsal:

```bash
pixi install --locked -e release
pixi run --locked -e release release-check
pixi run --locked -e release release-build
pixi run --locked -e release python scripts/release.py conda --sdist dist/gmri2fem-0.3.1.tar.gz
```

Use a clean artifact directory for each rehearsal. Release commands do not
publish artifacts; Conda packages are assembled into a local channel under
`dist/`. Wheel and Conda builds use the same staged source distribution.

The legacy `environment.yml` has been retired. Use Pixi for this library and
Gonzo's own environment for the full pipeline. `uv.lock` is a historical
Python-only snapshot, not the supported native environment or release
toolchain.

See the [release guide](docs/releases.md) for publication gates, companion-only
updates, artifact retries, and dependency policy.
