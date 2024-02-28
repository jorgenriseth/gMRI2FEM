import argparse
import logging
import os
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("--subjectid", required=True)
parser.add_argument("--t1file", required=True)
parser.add_argument("--t2file", default=None)
parser.add_argument("--outputdir", required=True)
parser.add_argument("--subjectsdir", default=f"{os.environ.get('FREESURFER')}/subjects")
parser.add_argument("-parallel", action="store_true")
args = parser.parse_args()

recon_all_cmd = (
    "recon-all"
    + f" -sd {args.subjectsdir}"
    + f" -s '{args.subjectid}'"
    + f" -i '{args.t1file}'"
    + (args.t2file is not None) * f" -T2 '{args.t2file}' -T2pial"
    + args.parallel * " -parallel"
    + " -all"
)
logger.info(f"Running cmd: '{recon_all_cmd}'")
subprocess.run(recon_all_cmd, shell=True).check_returncode()
shutil.move(Path(args.subjectsdir) / args.subjectid, Path(args.outputdir))
