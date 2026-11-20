import shutil

import click


def require_commands(*commands: str) -> None:
    missing = [command for command in commands if shutil.which(command) is None]
    if missing:
        raise click.ClickException(
            f"Missing external commands: {', '.join(missing)}. "
            "Install the gmri2fem-fsl Conda companion for FSL tools; "
            "install Greedy separately for reslicing. See the installation "
            "instructions at https://github.com/jorgenriseth/gMRI2FEM#installation."
        )
