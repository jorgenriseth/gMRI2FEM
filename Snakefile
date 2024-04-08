import os
import numpy as np
from pathlib import Path
from gmri2fem.filters import is_datetime_nii


shell.executable("/bin/bash")

configfile: "snakeconfig.yaml"

# To enable local scratch disks on clusters.
if workflow.run_local:
    workflow._shadow_prefix = os.environ.get("LOCALSCRATCH")

SESSIONS=[f"ses-{i+1:02d}" for i in range(config["num_sessions"])]
print(SESSIONS)
rule all:
  input:
    T1maps_LL = expand(
        "data/mri_dataset/derivatives/{subject}/{session}/anat/{subject}_{session}_T1map_LL_auto.nii.gz",
        subject=config["subjects"],
        session=SESSIONS
    )
    # expand(
    #     "data/mri_processed_data/modeling/resolution{res}/data.hdf",
    #     subject=config["subjects"],
    #     res=config["resolution"],
    #     session=SESSIONS,
    # )
    # expand(
    #     "data/mri_processed_data/{subject}/concentrations/{subject}_{session}_concentration_T1w.mgz",
    #     subject=config["subjects"],
    #     res=config["resolution"],
    #     session=SESSIONS,
    # )

rule T1map_estimation_from_LL:
  input:
    "data/mri_dataset/{subject}/{session}/anat"
  output:
    "data/mri_dataset/derivatives/{subject}/{session}/anat/{subject}_{session}_T1map_LL_auto.nii.gz"
  shell:
    "python gmri2fem/mriprocessing/looklocker_to_T1map.py"
    " --inputdir {input}"
    " --output {output}"
    " --threshold_number 3"
    # " --mask_quantile 0.1"


rule mri_convert_all:
  input:
    expand(
      "data/mri_processed_data/{subject}/conformed/{subject}_{session}_{sequence}_conformed.mgz",
      subject=config["subjects"],
      session=SESSIONS,
      sequence=["T1w", "T1map_LL_auto"]
    )


rule mri_convert:
    input:
        lambda wc: (
            "data/mri_dataset/derivatives/{subject}/{session}/anat/{subject}_{session}_{sequence}.nii.gz"
            if any([(x in wc.sequence) for x in ["T1map_LL_auto", "T1map_LookLocker"]]) else (
              "data/mri_dataset/{subject}/{session}/dwi/{subject}_{session}_{sequence}.nii.gz"
              if ("DTI" in wc.sequence) else
              "data/mri_dataset/{subject}/{session}/anat/{subject}_{session}_{sequence}.nii.gz"
            )
          )
    output:
        "data/mri_processed_data/{subject}/conformed/{subject}_{session,[A-Za-z0-9\-]+}_{sequence}_conformed.mgz"
    shell:
        "mri_convert --conform -odt float {input} {output}"


rule template_creation:
    input:
        expand(
            "data/mri_processed_data/{{subject}}/conformed/{{subject}}_{session}_T1w_conformed.mgz",
            session=SESSIONS
        )
    output:
        "data/mri_processed_data/{subject}/{subject}_T1w_template.mgz"
    shell:
        "mri_robust_template"
        " --mov {input}"
        " --template {output}"
        " --satit"
        " --inittp 1"
        " --fixtp"
        " --average 0"
        " --iscale"


ruleorder: registration_T1w > registration_bimodal

rule registration_T1w:
    input:
        image="data/mri_processed_data/{subject}/conformed/{subject}_{session}_T1w_conformed.mgz",
        template="data/mri_processed_data/{subject}/{subject}_T1w_template.mgz"
    output:
        image = "data/mri_processed_data/{subject}/registered/{subject}_{session}_T1w_registered.mgz",
        lta = "data/mri_processed_data/{subject}/registered/{subject}_{session}_T1w_registered.lta"
    shell:
        "mri_robust_register"
        " --mov {input.image}"
        " --dst {input.template}"
        " --lta {output.lta}"
        " --mapmov {output.image}"
        " --satit"
        " --iscale"

rule registration_bimodal:
    input:
        image="data/mri_processed_data/{subject}/conformed/{subject}_{session}_{sequence}_conformed.mgz",
        template="data/mri_processed_data/{subject}/{subject}_T1w_template.mgz"
    output:
        image = "data/mri_processed_data/{subject}/registered/{subject}_{session}_{sequence}_registered.mgz",
        lta = "data/mri_processed_data/{subject}/registered/{subject}_{session}_{sequence}_registered.lta"
    resources:
        time="05:00:00",
    shadow:
        config["shadow"]
    shell:
        "mri_robust_register"
        " --mov {input.image}"
        " --dst {input.template}"
        " --lta {output.lta}"
        " --mapmov {output.image}"
        " --cost NMI"


rule concentration_estimate:
    input:
        image="data/mri_processed_data/{subject}/registered/{subject}_{session}_T1map_LL_auto_registered.mgz",
        reference="data/mri_processed_data/{subject}/registered/{subject}_ses-01_T1map_LL_auto_registered.mgz",
    output:
        "data/mri_processed_data/{subject}/concentrations/{subject}_{session}_concentration_LL.mgz"
    shell:
        "python gmri2fem/mriprocessing/estimate_concentration.py"
        " --input {input.image}"
        " --reference {input.reference}"
        " --output {output}"
        " --r1 0.0032"


rule grow_refroi:
    input:
        seed="data/mri_processed_data/{subject}/{subject}_refroi_seed.mgz",
        references=expand(
            "data/mri_processed_data/{{subject}}/registered/{{subject}}_{session}_T1w_registered.mgz",
            session=SESSIONS
        )
    output:
        "data/mri_processed_data/{subject}/{subject}_refroi.mgz",
    shell:
        "python gmri2fem/mriprocessing/grow_refroi.py"
        " --seed {input.seed}"
        " --references {input.references}"
        " --output {output}"


rule normalize_T1w:
    input:
        image="data/mri_processed_data/{subject}/registered/{subject}_{session}_T1w_registered.mgz",
        refroi="data/mri_processed_data/{subject}/{subject}_refroi.mgz",
    output:
        "data/mri_processed_data/{subject}/T1w_normalized/{subject}_{session}_T1w_normalized.mgz"
    shell:
        "python gmri2fem/mriprocessing/normalize_images.py"
        " --image {input.image}"
        " --refroi {input.refroi}"
        " --output {output}"


rule T1map_literature:
    input:
        "data/mri_processed_data/freesurfer/{subject}/mri/aseg.mgz"
    output:
        t1map_synth = "data/mri_processed_data/{subject}/t1map_literature.mgz",
        brainmask = "data/mri_processed_data/{subject}/brainmask.mgz"
    shell:
        "python gmri2fem/mriprocessing/t1maps.py"
        " --aseg {input}"
        " --t1map {output.t1map_synth}"
        " --mask {output.brainmask}"

rule T1maps_T1w_estimated:
    input:
        T1w="data/mri_processed_data/{subject}/T1w_normalized/{subject}_{session}_T1w_normalized.mgz",
        T1w0="data/mri_processed_data/{subject}/T1w_normalized/{subject}_ses-01_T1w_normalized.mgz",
        T1map0="data/mri_processed_data/{subject}/t1map_literature.mgz"
    output:
        "data/mri_processed_data/{subject}/T1map_T1w/{subject}_{session}_T1map_T1w.mgz"
    shell:
        "python gmri2fem/mriprocessing/estimate_t1maps.py"
        " --inputT1w {input.T1w}"
        " --referenceT1w {input.T1w0}"
        " --referenceT1map {input.T1map0}"
        " --output {output}"


rule concentration_estimate_T1w: 
    input:
        image="data/mri_processed_data/{subject}/T1map_T1w/{subject}_{session}_T1map_T1w.mgz",
        reference="data/mri_processed_data/{subject}/T1map_T1w/{subject}_ses-01_T1map_T1w.mgz",
        mask="data/mri_processed_data/{subject}/brainmask.mgz",
    output:
        "data/mri_processed_data/{subject}/concentrations/{subject}_{session}_concentration_T1w.mgz"
    shell:
        "python gmri2fem/mriprocessing/estimate_concentration.py"
        " --input {input.image}"
        " --reference {input.reference}"
        " --output {output}"


rule recon_all_base:
    input:
        t1="data/mri_processed_data/{subject}/conformed/{subject}_ses-01_T1w_conformed.mgz",
    output:
        protected(directory("data/mri_processed_data/freesurfer/{subject}/")),
    resources:
        time="24:00:00",
    shadow:
        config["shadow"]
    shell:
        "recon-all"
        " -sd 'data/mri_processed_data/freesurfer'"
        " -s {wildcards.subject}"
        " -i {input.t1}"
        " -all"


rule recon_all_FLAIR:
    input:
        t1="data/mri_processed_data/{subject}/conformed/{subject}_ses-01_T1w_conformed.mgz",
        FLAIR="data/mri_processed_data/{subject}/conformed/{subject}_ses-01_FLAIR_conformed.mgz",
    output:
        protected(directory("freesurfer/{subject}_FLAIR/")),
    resources:
        time="48:00:00",
    shadow:
        config["shadow"]
    shell:
        "recon-all"
        " -sd 'freesurfer'"
        " -s {wildcards.subject}_FLAIR"
        " -i {input.t1}"
        " -FLAIR {input.FLAIR}"
        " -FLAIRpial"
        " -all"


rule recon_all_T2w:
    input:
        t1="data/mri_processed_data/{subject}/conformed/{subject}_ses-01_T1w_conformed.mgz",
        t2="data/mri_processed_data/{subject}/conformed/{subject}_ses-01_T2w_conformed.mgz",
    output:
        protected(directory("freesurfer/{subject}_T2w/")),
    resources:
        time="96:00:00",
    shadow:
        config["shadow"]
    shell:
        "recon-all"
        " -sd 'freesurfer'"
        " -s {wildcards.subject}_T2w"
        " -i {input.t1}"
        " -T2 {input.t2}"
        " -T2pial"
        " -all"


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

rule dilated_mask:
    input:
        "data/mri_processed_data/freesurfer/{subjectid}/mri/wmparc.mgz"
    output: 
        "data/mri_processed_data/{subjectid}/mask_dilated.mgz"
    shell:
        "mri_binarize --i {input} --gm --dilate 2 --o {output}"

rule clean_dti:
    input:
        dtifile = "data/mri_processed_data/freesurfer/{subjectid}/dti/tensor.nii.gz",
        maskfile = "data/mri_processed_data/{subjectid}/mask_dilated.mgz"
    output:
        dtifile = "data/mri_processed_data/{subjectid}/dti/{subjectid}_dti-tensor_clean.mgz",
    shell:
        "python gmri2fem/dti/clean_dti_data.py"
        " --dti {input.dtifile}"
        " --mask {input.maskfile}"
        " --out {output}"

rule dti2hdf:
    input:
        meshfile = "data/mri_processed_data/{subjectid}/modeling/resolution32/mesh.hdf",
        dtifile= "data/mri_processed_data/{subjectid}/dti/{subjectid}_dti-tensor_clean.mgz"
    output:
        "data/mri_processed_data/{subjectid}/modeling/resolution{res}/dti.hdf"
    shell:
        "python gmri2fem/dti/dti_data_to_mesh.py"
        " --dti {input.dtifile}"
        " --mesh {input.meshfile}"
        " --out {output}"



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
