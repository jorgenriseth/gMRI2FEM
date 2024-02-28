import os
import numpy as np
from pathlib import Path
from gmri2fem.filters import is_datetime_nii


shell.executable("/bin/bash")

configfile: "snakeconfig.yaml"

# To enable local scratch disks on clusters.
if workflow.run_local:
    workflow._shadow_prefix = os.environ.get("LOCALSCRATCH")


SESSIONS=[f"ses-{i+1:02d}" for i in range(5)]
rule all:
  input:
    expand(
      "data/mri_processed_data/sub-01/conformed/sub-01_{session}_{sequence}_conformed.mgz",
      session=SESSIONS, sequence="T1w")

rule mri_convert:
    input:
        lambda wc: (
            "data/mri_dataset/derivatives/{subject}/{session}/anat/{subject}_{session}_{sequence}.nii.gz"
            if ("T1map_LookLocker" in wc.sequence) else
            "data/mri_dataset/{subject}/{session}/dwi/{subject}_{session}_{sequence}.nii.gz"
            if ("DTI" in wc.sequence) else
            "data/mri_dataset/{subject}/{session}/anat/{subject}_{session}_{sequence}.nii.gz"
        )
    output:
        "data/mri_processed_data/{subject}/conformed/{subject}_{session}_{sequence}_conformed.mgz"
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
        image = "data/mri_processed_data/{subject}/{sequence}_registered/{subject}_{session}_T1w_registered.mgz",
        lta = "data/mri_processed_data/{subject}/{sequence}_registered/{subject}_{session}_T1w_registered.lta"
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
        image = "data/mri_processed_data/{subject}/{sequence}_registered/{subject}_{session}_{sequence}_registered.mgz",
        lta = "data/mri_processed_data/{subject}/{sequence}_registered/{subject}_{session}_{sequence}_registered.lta"
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
        image="data/mri_processed_data/{subject}/T1map_LookLocker_registered/{subject}_{session}_T1map_LookLocker_registered.mgz",
        reference="data/mri_processed_data/{subject}/T1map_LookLocker_registered/{subject}_ses-01_T1map_LookLocker_registered.mgz",
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
            "data/mri_processed_data/{{subject}}/T1w_registered/{{subject}}_{session}_T1w_registered.mgz",
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
        image="data/mri_processed_data/{subject}/T1w_registered/{subject}_{session}_T1w_registered.mgz",
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
        "freesurfer/{subject}/mri/aseg.mgz"
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
        reference="data/mri_processed_data/{subject}/T1w_normalized/{subject}_{session}_T1w_normalized.mgz",
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
        protected(directory("freesurfer/{subject}/")),
    resources:
        time="24:00:00",
    shadow:
        config["shadow"]
    shell:
        "recon-all"
        " -sd 'freesurfer'"
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

#rule extract_times_workflow:
#    input:
#        injection="data/injection_times.json",
#        t1dir="data/{subject}/RESAMPLED",
#    output:
#        "data/{subject}/timestamps.txt",
#    shell:
#        "python gmri2fem/mriprocessing/extract_timestamps.py"
#        " --subject {wildcards.subject}"
#        " --injection_times {input.injection}"
#        " --t1dir {input.t1dir}"
#        " --output {output}"
#
#
#rule region_statistics_workflow:
#    input:
#        timestamps="data/{subject}/timestamps.txt",
#        asegfile="data/{subject}/freesurfer/mri/aseg.mgz",
#        concentrationdir="data/{subject}/CONCENTRATIONS",
#    output:
#        "data/{subject}/STATISTICS/lut-regions-stats.csv",
#    shell:
#        "python gmri2fem/analysis/region_quantities.py "
#        " --asegfile {input.asegfile}"
#        " --timestamps {input.timestamps}"
#        " --concentrationdir {input.concentrationdir}"
#        " --lutfile data/freesurfer_lut.json"
#        " --output {output}"
#
#
#rule region_statistics:
#    input:
#        expand(
#            "data/{subject}/STATISTICS/lut-regions-stats.csv",
#            subject=config["subjects"],
#        ),
#        T1files=[
#            Path(f"data/{subject}/RESAMPLED/{file.with_suffix('.mgz').name}")
#            for subject, files in SUBJECT_FILES.items()
#            for file in files["T1"]
#        ],
#        T2file=[
#            Path(f"data/{subject}/RESAMPLED/T2.mgz") for subject in config["subjects"]
#        ],
#
#
#rule concentration_estimate_workflow:
#    input:
#        image="data/{subjectid}/T1MAPS/{image}",
#        reference=lambda x: f"data/{x.subjectid}/T1MAPS/{T1_FILES[x.subjectid][0].stem}.mgz",
#        mask="data/{subjectid}/brainmask.mgz",
#    output:
#        "data/{subjectid}/CONCENTRATIONS/{image}"
#    shell:
#        "python gmri2fem/mriprocessing/estimate_concentration.py"
#        " --input {input.image}"
#        " --reference {input.reference}"
#        " --output {output}"
#        # " --mask {input.mask}"
#
#rule concentration_estimate:
#    input:
#        T1files=[
#            Path(f"data/{subject}/CONCENTRATIONS/{file.with_suffix('.mgz').name}")
#            for subject, files in SUBJECT_FILES.items()
#            for file in files["T1"]
#        ],
#
#
#rule recon_all_workflow:
#    input:
#        t1=lambda x: sorted(Path(f"data/{x.subjectid}/RESAMPLED").iterdir())[0],
#        t2="data/{subjectid}/RESAMPLED/T2.mgz",
#        script="mri_recon.py",
#    output:
#        protected(directory("data/freesurfer/{subjectid}/")),
#    threads: config["recon_threads"]
#    resources:
#        time="96:00:00",
#    params:
#        subjectsdir=config["FS_subjectsdir"],
#    shadow:
#        config["shadow"]
#    shell:
#        "python {input.script}"
#        " --subjectid {wildcards.subjectid}"
#        " --t1file {input.t1}"
#        " --outputdir {output}"
#        " --subjectsdir {params.subjectsdir}"
#        # " -parallel"
#        # " --t2file {input.t2}"
#
#
#rule symlink_freesurfer_workflow:
#    input:
#        "data/freesurfer/{subjectid}",
#    output:
#        directory("data/{subjectid}/freesurfer"),
#    shell:
#        "ln -s '../../{input}' '{output}'"
#
#
#rule mri_recon_all:
#    input:
#        expand("data/{subjectid}/freesurfer", subjectid=config["subjects"]),
#
#
#####################################
## Converting data to FEniCS-formats.
#####################################
#
#rule mri2fenics:
#    input:
#        meshfile="data/{subjectid}/MODELING/resolution{res}/mesh.hdf",
#        timestamps="data/{subjectid}/timestamps.txt",
#        concentration_files= lambda wc: [
#            f"data/{wc.subjectid}/CONCENTRATIONS/{file.with_suffix('.mgz').name}"
#              for file in T1_FILES[wc.subjectid]
#        ],
#    output:
#        hdf="data/{subjectid}/MODELING/resolution{res}/data.hdf",
#        visual="data/{subjectid}/MODELING/resolution{res}/visual/data.xdmf",
#    shell:
#        "python gmri2fem/mriprocessing/mri2fenics.py"
#        " --mris {input.concentration_files}"
#        " --mesh_hdf {input.meshfile}"
#        " --output_hdf {output.hdf}"
#        " --timestamps {input.timestamps}"
#
#rule mri2fenics_all:
#    input:
#        expand(
#            "data/{subjectid}/MODELING/resolution{res}/data.hdf",
#            subjectid=config["subjects"],
#            res=config["resolution"],
#        ),
#
#
#rule surface_stl_conversion:
#    input:
#        "data/{subjectid}/freesurfer/surf/{filename}",
#    output:
#        "data/{subjectid}/MODELING/surfaces/{filename}.stl",
#    shell:
#        "mris_convert --to-scanner {input} {output}"
#
#
#rule ventricle_extraction:
#    input:
#        "data/{subjectid}/freesurfer/mri/wmparc.mgz",
#    output:
#        fs="data/{subjectid}/MODELING/surfaces/ventricles",
#        stl="data/{subjectid}/MODELING/surfaces/ventricles.stl",
#    shell:
#        "bash scripts/extract-ventricles.sh {input} {output.fs}"
#        " && mris_convert --to-scanner {output.fs} {output.stl}"
#
#
#SURFACES = ["lh.pial", "rh.pial", "lh.white", "rh.white", "ventricles"]
#
#
#rule create_surfaces:
#    input:
#        expand(
#            "data/{subjectid}/MODELING/surfaces/{filename}.stl",
#            filename=SURFACES,
#            subjectid=config["subjects"],
#        ),
#
#
#rule create_mesh:
#    input:
#        expand(
#            "data/{{subjectid}}/MODELING/surfaces/{filename}.stl",
#            filename=SURFACES,
#        ),
#    output:
#        "data/{subjectid}/MODELING/resolution{res}/mesh.hdf",
#    resources:
#        time="00:30:00",
#    shadow:
#        config["shadow"]
#    shell:
#        "python gmri2fem/mriprocessing/mesh_generation.py"
#        " --surfaces {input}"
#        " --output {output}"
#        " --resolution {wildcards.res}"
#
#
############################
## Simulation postprocessing
############################
#rule fenics2mri_workflow:
#    input:
#        referenceimage="data/{subjectid}/freesurfer/mri/T1.mgz",
#        simulationfile="data/{subjectid}/MODELING/resolution{res}/{funcname}.hdf",
#        timestampfile="data/{subjectid}/timestamps.txt",
#    output:
#        "data/{subjectid}/MODELING/resolution{res}/MRIs/{funcname}_{idx}.nii.gz",
#    shell:
#        "python gmri2fem/mriprocessing/fenics2mri.py"
#        " --simulationfile {input.simulationfile}"
#        " --output {output}"
#        " --referenceimage {input.referenceimage}"
#        " --timestamps {input.timestampfile}"
#        " --timeidx {wildcards.idx}"
#        " --functionname 'total_concentration'"
#
#rule solute_quantification_workflow:
#    input:
#        "data/{subjectid}/MODELING/resolution{res}/{funcname}.hdf",
#    output:
#        "data/{subjectid}/MODELING/resolution{res}/quantities_{funcname}.csv",
#    threads: config["sim_threads"]
#    shadow:
#        config["shadow"]
#    shell:
#        "OMP_NUM_THREADS=1 mpirun -n {threads}"
#        " python gmri2fem/analysis/solute_quantification.py"
#        " --input {input} --funcname total_concentration --output {output}"
#
#
#rule fenics2mri_transpose:
#    input:
#        expand(
#            "data/{subjectid}/MODELING/resolution{res}/MRIs/data_{idx}.nii.gz",
#            subjectid=config["subjects"],
#            res=config["resolution"],
#            idx=range(5)
#        )
#
