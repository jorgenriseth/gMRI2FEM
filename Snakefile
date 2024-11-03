import os
import numpy as np
from pathlib import Path

singularity: "singularity/gmri2fem.sif"
shell.executable("bash")

configfile: "snakeconfig.yaml"

if DeploymentMethod.APPTAINER in workflow.deployment_settings.deployment_method:
  shell.prefix(
    "set -eo pipefail; "
    + "source /opt/conda/etc/profile.d/conda.sh && "
    + "conda activate $CONDA_ENV_NAME && "
    #    + "ls -l && "
  )

wildcard_constraints:
  session = "ses-\d{2}"


if "subjects" in config:
	SUBJECTS=config["subjects"]
else:
  SUBJECTS = [p.stem for p in Path("data/mri_dataset/sourcedata").glob("sub-*")]
  config["subjects"] = SUBJECTS

if "ignore_subjects" in config:
  for subject in config["ignore_subjects"]:
    if subject in SUBJECTS:
      SUBJECTS.remove(subject)
  config["subjects"] = SUBJECTS

SESSIONS = {
  subject: sorted([p.stem for p in Path(f"data/mri_dataset/sourcedata/{subject}").glob("ses-*")])
  for subject in SUBJECTS
}
config["sessions"] = SESSIONS


include: "workflows/T1maps"
include: "workflows/register"
#include: "workflows/recon-all"
include: "workflows/T1w_signal_intensities"
include: "workflows/concentration-estimate"
include: "workflows/statistics"
include: "workflows/mesh-generation"
include: "workflows/mri2fem"
include: "workflows/dti"
include: "workflows/snakefile-modeling"

