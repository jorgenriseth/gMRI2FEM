import os
import numpy as np
from pathlib import Path

singularity: "singularity/gmri2fem.sif"
shell.executable("bash")

configfile: "snakeconfig.yaml"

if workflow.use_singularity:
  shell.prefix(
    "set -eo pipefail; "
    + "source /opt/conda/etc/profile.d/conda.sh && "
    + "conda activate $CONDA_ENV_NAME && "
    #    + "ls -l && "
  )

# To enable local scratch disks on clusters.
# Not sure if works properly.
if workflow.run_local:
    workflow._shadow_prefix = os.environ.get("LOCALSCRATCH")

if workflow._shadow_prefix is not None:
  shadow_directive = "full"


wildcard_constraints:
  session = "ses-\d{2}"

SESSIONS=[f"ses-{i+1:02d}" for i in range(config["num_sessions"])]
rule all:
  input:
    T1maps_LL = expand(
        "data/mri_dataset/derivatives/{subject}/{session}/anat/{subject}_{session}_T1map_LL_auto.nii.gz",
        subject=config["subjects"],
        session=SESSIONS
    )



module preprocessing:
  snakefile: "data/mri_dataset/Snakefile"
  prefix: "data/mri_dataset"
  config: config

use rule * from preprocessing as preprocessing_*

include: "workflows_additional/register"
include: "workflows_additional/recon-all"
include: "workflows_additional/T1maps"
include: "workflows_additional/concentration-estimate"
include: "workflows_additional/statistics"
include: "workflows_additional/mesh-generation"
include: "workflows_additional/mri2fem"
