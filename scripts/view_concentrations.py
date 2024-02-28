#!/usr/bin/env python3
import subprocess
from pathlib import Path


def set_colormap(path: Path, visible: bool = False) -> str:
    return f"{path}:colormap=turbo:colorscale=0,1:visible={int(visible)}"


currentdir = Path(".")
filelist = sorted(currentdir.iterdir())

arguments = [set_colormap(filelist[0], True), *[set_colormap(p) for p in filelist[1:]]]
arguments.reverse()
cmd = f"freeview \
    {' '.join(arguments)}"
print(cmd)
subprocess.run(cmd, shell=True)
