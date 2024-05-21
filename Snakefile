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
  )

# To enable local scratch disks on clusters.
if workflow.run_local:
    workflow._shadow_prefix = os.environ.get("LOCALSCRATCH")

if workflow._shadow_prefix is not None:
  shadow_directive = "full"


SESSIONS=[f"ses-{i+1:02d}" for i in range(config["num_sessions"])]
print(SESSIONS)

wildcard_constraints:
  session = "ses-\d{2}"

rule all:
  input:
    T1maps_LL = expand(
        "data/mri_dataset/derivatives/{subject}/{session}/anat/{subject}_{session}_T1map_LL_auto.nii.gz",
        subject=config["subjects"],
        session=SESSIONS
    )


include: "data/mri_dataset/Snakefile"
include: "workflows_additional/register"
include: "workflows_additional/recon-all"
include: "workflows_additional/T1maps"
include: "workflows_additional/concentration-estimate"
include: "workflows_additional/statistics"


rule extract_concentration_times_T1w:
   input:
       timetable="data/mri_dataset/derivatives/timetable.csv"
   output:
       "data/mri_processed_data/{subject}/timestamps_T1w.txt",
   shell:
       "python gmri2fem/mriprocessing/extract_timestamps.py"
       " --timetable {input}"
       " --subject {wildcards.subject}"
       " --sequence_label T1w"
       " --output {output}"


rule extract_concentration_times_LL:
   input:
       timetable="data/mri_dataset/derivatives/timetable.csv"
   output:
       "data/mri_processed_data/{subject}/timestamps_LL.txt"
   shell:
       "python gmri2fem/mriprocessing/extract_timestamps.py"
       " --timetable {input}"
       " --subject {wildcards.subject}"
       " --sequence_label LookLocker"
       " --output {output}"

rule mri2fenics:
    input:
        meshfile="data/mri_processed_data/{subjectid}/modeling/resolution{res}/mesh.hdf",
        timestamps="data/mri_processed_data/{subjectid}/timestamps_LL.txt",
        concentrations=expand(
          "data/mri_processed_data/{{subjectid}}/concentrations/{{subjectid}}_{session}_concentration_LL.mgz",
          session=SESSIONS
        )
    output:
        hdf="data/mri_processed_data/{subjectid}/modeling/resolution{res}/data.hdf",
        visual="data/mri_processed_data/{subjectid}/modeling/resolution{res}/visual/data.xdmf",
    shell:
        "python gmri2fem/mriprocessing/mri2fenics.py"
        " --mris {input.concentrations}"
        " --mesh_hdf {input.meshfile}"
        " --output_hdf {output.hdf}"
        " --timestamps {input.timestamps}"



rule surface_stl_conversion:
    input:
        "data/mri_processed_data/freesurfer/{subjectid}/surf/{filename}",
    output:
        "data/mri_processed_data/{subjectid}/modeling/surfaces/{filename}.stl",
    shell:
        "mris_convert --to-scanner {input} {output}"


rule ventricle_extraction:
    input:
        "data/mri_processed_data/freesurfer/{subjectid}/mri/wmparc.mgz",
    output:
        fs="data/mri_processed_data/{subjectid}/modeling/surfaces/ventricles",
        stl="data/mri_processed_data/{subjectid}/modeling/surfaces/ventricles.stl",
    shell:
        "bash scripts/extract-ventricles.sh {input} {output.fs}"
        " && mris_convert --to-scanner {output.fs} {output.stl}"


SURFACES = ["lh.pial", "rh.pial", "lh.white", "rh.white", "ventricles"]
rule create_surfaces:
    input:
        expand(
            "data/mri_processed_data/{{subjectid}}/modeling/surfaces/{filename}.stl",
            filename=SURFACES,
        ),


rule create_mesh:
    input:
        expand(
            "data/mri_processed_data/{{subjectid}}/modeling/surfaces/{filename}.stl",
            filename=SURFACES,
        ),
    output:
        "data/mri_processed_data/{subjectid}/modeling/resolution{res}/mesh.hdf",
    resources:
        time="00:30:00",
    shadow:
        config["shadow"]
    shell:
        "python gmri2fem/mriprocessing/mesh_generation.py"
        " --surfaces {input}"
        " --output {output}"
        " --resolution {wildcards.res}"


rule region_statistics_workflow:
   input:
       timestamps="data/mri_processed_data/{subject}/timestamps_LL.txt",
       segfile="data/mri_processed_data/freesurfer/{subject}/mri/aseg.mgz",
       concentrationdir="data/mri_processed_data/{subject}/concentrations",
   output:
       "data/mri_processed_data/{subject}/analysis/lut-regions-stats.csv",
   shell:
       "python gmri2fem/analysis/region_quantities.py "
       " --asegfile {input.asegfile}"
       " --timestamps {input.timestamps}"
       " --concentrationdir {input.concentrationdir}"
       " --lutfile data/mri_processed_data/freesurfer_lut.json"
       " --output {output}"


rule fenics2mri_workflow:
   input:
       referenceimage="data/mri_processed_data/freesurfer/{subjectid}/mri/T1.mgz",
       simulationfile="data/mri_processed_data/{subject}/modelingn/resolution{res}/{funcname}.hdf",
       timestampfile="data/{subjectid}/timestamps_LL.txt",
   output:
       "data/{subjectid}/MODELING/resolution{res}/MRIs/{funcname}_{idx}.nii.gz",
   shell:
       "python gmri2fem/mriprocessing/fenics2mri.py"
       " --simulationfile {input.simulationfile}"
       " --output {output}"
       " --referenceimage {input.referenceimage}"
       " --timestamps {input.timestampfile}"
       " --timeidx {wildcards.idx}"
       " --functionname 'total_concentration'"


rule fenics2mri_transpose:
   input:
       expand(
           "data/{subjectid}/MODELING/resolution{res}/MRIs/data_{idx}.nii.gz",
           subjectid=config["subjects"],
           res=config["resolution"],
           idx=range(5)
       )
