# Building and releasing

## Authoritative inputs

- Python metadata: [pyproject.toml](../pyproject.toml).
- Conda-only translations for PyVista extras and tested version overrides:
  [dependencies.toml](../conda-recipe/dependencies.toml).
- FSL compatibility policy and independent companion version/build:
  [fsl.toml](../conda-recipe/fsl.toml).
- Exact contributor/release-tool environments: [pixi.lock](../pixi.lock).

The Conda recipe derives ordinary Python requirements directly from
`pyproject.toml`; only dependencies needing Conda-specific expansion or version
overrides are listed separately in `dependencies.toml`. The release check
requires every override to correspond to a declared Python dependency and to
have a documented reason. The recipe deliberately supports PyVista 0.44 and
expands its `all` extra explicitly, matching Gonzo's compatibility cap. Newer
PyVista extras differ: review that override and test the pipeline before
widening it.
VTK retains Gonzo's tested cap below 9.4: PyVista 0.44 imports
`vtkCapsuleSource`, which newer VTK versions removed.
The Conda host and native test interpreter are deliberately Python 3.11;
runtime metadata retains the project's broader Python range. The wheel is
tested separately on both supported Python versions.
The FEniCS mesh/function HDF5 I/O formerly provided by panta-rhei lives in
the bundled `fenics_io` package; panta-rhei (and its Gmsh requirement, which
previously failed the strict `pip check` gate) is no longer a dependency.
Keep `pip check` strict: do not patch installed metadata or add checker
exceptions. Local build artifacts are rehearsal outputs, not approved releases.

FSL utilities `fslroi`, `fslmerge`, and `fslmaths` are supplied by
`fsl-avwutils`. Topup additionally requires its `b02b0.cnf` configuration.
The companion does not install Greedy or license-restricted FreeSurfer.

## Rehearsal (no publication)

Use Pixi 0.83.0 and a fresh checkout/artifact directory:

```bash
pixi run --locked -e release python scripts/release.py check
pixi run --locked -e release python -m pytest -q tests/test_release.py
pixi run --locked -e release python scripts/release.py build
pixi run --locked -e release python scripts/release.py conda --sdist dist/gmri2fem-0.3.1.tar.gz
```

Substitute the actual project version for the archive filename. The wheel is
built from the sdist, not from a separate checkout. The Conda recipe is copied
beside that same extracted source. The companion resolves against the local
channel, so the new library need not already exist on PyPI or Anaconda.

Builds preserve checksums and commit provenance in `dist`. Install the wheel
in a fresh environment outside the checkout and run `test_cli.py`,
`test_concentration.py`, and `test_reslice_dti.py` from the sdist. Conda builds
use `--no-test` for packaging, then the orchestrator uses Pixi to install the
exact-version artifacts from the local channel and execute the native/CLI and
FSL utility tests. This avoids Conda-build's slow runtime solver without
skipping validation. Any failed install or test fails the build command and
blocks publication. The recipes retain tests for standalone Conda-build use.
Test manifests/locks are preserved in `dist/test-environments` for provenance;
their local-channel paths are not portable consumer environments.
Real-data scientific comparisons remain a separate promotion gate.

## Publication

The release workflow replaces both old independent upload workflows. PRs,
main pushes, and default manual runs build/test only. Publication requires a
matching tag and either a published GitHub release or explicit manual opt-in.
Before enabling it:

1. Select an unused upstream version and matching `vVERSION` tag; do not
   overwrite existing PyPI releases.
2. Configure PyPI trusted publishing for the `pypi` GitHub environment and
   require environment approval.
3. Configure `ANACONDA_TOKEN` in the protected `anaconda` environment.
4. Approve only after wheel tests on both Python versions, Conda/native tests,
   companion tests, and downstream Gonzo checks pass.

Publication jobs download the exact tested artifacts, verify checksums, and
upload without `--force` or `--skip-existing`. Conda uploads the library before
the companion. PyPI/Anaconda do not provide a shared transaction.

For a partial upload failure, retain/download the workflow artifacts. Compare
the existing registry artifact's checksum before deciding that its upload is
complete, then manually upload only missing artifacts after approval. Do not
rebuild under an already published filename. If the artifact differs, use a
new Python version or Conda build number as appropriate.

## Companion-only updates

Change the FSL policy/version (or its build number for packaging corrections),
run the same rehearsal with `--companion-only` on the `conda` command, and
verify that the compatible library already exists in the configured channel.
Use tag `fsl-vVERSION` and a release/manual run with `companion_only` enabled.
Python artifacts may be built for source validation but are not published;
the Python project's version does not need to change.

Do not split more feature metapackages until there is an actual installation
need. Future companions should each depend directly on compatible gmri2fem
ranges, not require lockstep releases. Conda checks their constraint
intersection; tests must still establish combined behavioral compatibility.
