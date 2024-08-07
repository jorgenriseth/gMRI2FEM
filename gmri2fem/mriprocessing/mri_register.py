import logging
import subprocess
from pathlib import Path
from typing import Iterator

from gmri2fem.filters import is_T1_mgz

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def as_str(path_iter: Iterator[Path]) -> str:
    return " ".join(map(lambda x: f"'{x}'", path_iter))


def patient_create_template(
    input_dir: Path, registered_dir: Path, template: Path, lta_dir: Path
) -> Path:
    volumes = sorted(filter(is_T1_mgz, input_dir.iterdir()))
    template.parent.mkdir(exist_ok=True)
    lta_dir.mkdir(exist_ok=True)
    registered_dir.mkdir(exist_ok=True)
    template_command = (
        f"mri_robust_template"
        + f" --mov {as_str(volumes)}"
        + f" --template '{template}'"
        + " --satit"
        + f" --lta {as_str(map(lambda x: lta_dir / x.with_suffix('.lta').name, volumes))}"
        + f" --mapmov {as_str(map(lambda x: registered_dir / x.name, volumes))}"
        + " --average 1"
        + " --inittp 1"
        + " --fixtp"
        + " --iscale"
        + " --maxit 10"
        + " --highit 10"
    )
    print(template_command)
    subprocess.run(template_command, shell=True).check_returncode()
    return template


def patient_register(
    input_dir: Path, registered_dir: Path, template: Path, lta_dir: Path
) -> Path:
    lta_dir.mkdir(exist_ok=True)
    volumes = sorted(filter(is_T1_mgz, input_dir.iterdir()))
    for volume in volumes:
        register_command = (
            f"mri_robust_register"
            + f" --mov '{volume}'"
            + f" --dst '{template}'"
            + f" --mapmov '{registered_dir / volume.name}'"
            + f" --lta '{lta_dir / volume.with_suffix('.lta').name}'"
            + " --iscale"
            + " --satit "
            + " --maxit 10"
        )
        logger.info(register_command)
        subprocess.run(register_command, shell=True).check_returncode()
    return registered_dir


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", required=True, help="Directory to T1-weighted images to register."
    )
    parser.add_argument(
        "--output", required=True, help="Output directory for registered images."
    )
    parser.add_argument("--ltadir", required=True, help="LTA-directory storage path")
    args = parser.parse_args()

    indir = Path(args.input)
    template = Path(args.output) / "template.mgz"
    ltadir = Path(args.ltadir)
    outdir = Path(args.output)
    patient_create_template(indir, outdir, template, ltadir)
    # patient_register(indir, outdir, template, ltadir)
