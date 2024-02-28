import os
import numpy as np
from pathlib import Path
from gmri2fem.filters import is_datetime_nii


shell.executable("/bin/bash")

configfile: "snakeconfig.yaml"

# To enable local scratch disks on clusters.
if workflow.run_local:
    workflow._shadow_prefix = os.environ.get("LOCALSCRATCH")

SUBJECT_FILES = {
    subject: {
        "T1": sorted(filter(is_datetime_nii, Path(f"DATA/{subject}/T1/").iterdir())),
        "T2": next(filter(is_datetime_nii, Path(f"DATA/{subject}/T2").iterdir())),
    }
    for subject in config["subjects"]
}
T1_FILES = {
    subject: sorted(filter(is_datetime_nii, Path(f"DATA/{subject}/T1/").iterdir()))
    for subject in config["subjects"]
}

T1_FILESTEMS = {
    subject: file.stem 
    for subject in config["subjects"]
    for file in sorted(filter(is_datetime_nii, Path(f"DATA/{subject}/T1/").iterdir()))
}

rule mri_convert:
    input:
        T1files=[
            Path(f"DATA/{subject}/RESAMPLED/{file.with_suffix('.mgz').name}")
            for subject, files in SUBJECT_FILES.items()
            for file in files["T1"]
        ],
        T2file=[
            Path(f"DATA/{subject}/RESAMPLED/T2.mgz") for subject in config["subjects"]
        ],
    shadow:
        config["shadow"]


rule T1convert:
    input:
        "{subjectdata}/T1/{filename}.nii",
    output:
        "{subjectdata}/RESAMPLED/{filename}.mgz",
    shadow:
        config["shadow"]
    shell:
        "mri_convert --conform -odt float {input} {output}"


rule T2convert:
    input:
        lambda x: Path(f"{x.subjectdata}/T2").glob("*.nii"),
    output:
        "{subjectdata}/RESAMPLED/T2.mgz",
    shadow:
        config["shadow"]
    shell:
        "mri_convert --conform -odt float {input} {output}"


rule register:
    input:
        script="gmri2fem/mriprocessing/mri_register.py",
        data="{subject_data}/RESAMPLED/",
    output:
        reg=directory("{subject_data}/REGISTERED/"),
        lta=directory("{subject_data}/LTA/"),
    resources:
        time="05:00:00",
    shadow:
        config["shadow"]
    shell:
        "python {input.script} --input {input.data}  --output {output.reg} --ltadir {output.lta}"


rule register_all:
    input:
        expand("DATA/{subjectid}/REGISTERED/", subjectid=config["subjects"]),


rule t1map_reference:
    input:
        "DATA/{subjectid}/freesurfer/mri/aseg.mgz",
    output:
        t1="DATA/{subjectid}/t1map.mgz",
        mask="DATA/{subjectid}/brainmask.mgz",
    shell:
        "python gmri2fem/mriprocessing/t1maps.py"
        " --aseg {input} --t1map {output.t1} --mask {output.mask}"


rule t1map_estimation_workflow:
    input:
        t1w="DATA/{subjectid}/NORMALIZED/{t1wimage}",
        t1w0=lambda wc: f"DATA/{wc.subjectid}/NORMALIZED/{T1_FILES[wc.subjectid][0].stem}.mgz",
        t10="DATA/{subjectid}/t1map.mgz",
    output:
        "DATA/{subjectid}/T1MAPS/{t1wimage}"
    shell:
        "python gmri2fem/mriprocessing/estimate_t1maps.py"
        " --inputT1w {input.t1w}"
        " --referenceT1w {input.t1w0}"
        " --referenceT1map {input.t10}"
        " --output {output}"

rule t1map_estimation:
    input:
        T1files=[
            Path(f"DATA/{subject}/T1MAPS/{file.with_suffix('.mgz').name}")
            for subject, files in SUBJECT_FILES.items()
            for file in files["T1"]
        ]

rule mri_normalize_workflow:
    input:
        images="DATA/{subjectid}/REGISTERED",
        refroi="DATA/{subjectid}/refroi.mgz",
    output:
        directory("DATA/{subjectid}/NORMALIZED"),
    shell:
        "python gmri2fem/mriprocessing/normalize_images.py {wildcards.subjectid}"


rule mri_normalize:
    input:
        expand("DATA/{subjectid}/NORMALIZED", subjectid=config["subjects"]),


rule extract_times_workflow:
    input:
        injection="DATA/injection_times.json",
        t1dir="DATA/{subject}/RESAMPLED",
    output:
        "DATA/{subject}/timestamps.txt",
    shell:
        "python gmri2fem/mriprocessing/extract_timestamps.py"
        " --subject {wildcards.subject}"
        " --injection_times {input.injection}"
        " --t1dir {input.t1dir}"
        " --output {output}"


rule extract_times:
    input:
        expand("DATA/{subject}/timestamps.txt", subject=config["subjects"]),


rule region_statistics_workflow:
    input:
        timestamps="DATA/{subject}/timestamps.txt",
        asegfile="DATA/{subject}/freesurfer/mri/aseg.mgz",
        concentrationdir="DATA/{subject}/CONCENTRATIONS",
    output:
        "DATA/{subject}/STATISTICS/lut-regions-stats.csv",
    shell:
        "python gmri2fem/analysis/region_quantities.py "
        " --asegfile {input.asegfile}"
        " --timestamps {input.timestamps}"
        " --concentrationdir {input.concentrationdir}"
        " --lutfile DATA/freesurfer_lut.json"
        " --output {output}"


rule region_statistics:
    input:
        expand(
            "DATA/{subject}/STATISTICS/lut-regions-stats.csv",
            subject=config["subjects"],
        ),
        T1files=[
            Path(f"DATA/{subject}/RESAMPLED/{file.with_suffix('.mgz').name}")
            for subject, files in SUBJECT_FILES.items()
            for file in files["T1"]
        ],
        T2file=[
            Path(f"DATA/{subject}/RESAMPLED/T2.mgz") for subject in config["subjects"]
        ],


rule concentration_estimate_workflow:
    input:
        image="DATA/{subjectid}/T1MAPS/{image}",
        reference=lambda x: f"DATA/{x.subjectid}/T1MAPS/{T1_FILES[x.subjectid][0].stem}.mgz",
        mask="DATA/{subjectid}/brainmask.mgz",
    output:
        "DATA/{subjectid}/CONCENTRATIONS/{image}"
    shell:
        "python gmri2fem/mriprocessing/estimate_concentration.py"
        " --input {input.image}"
        " --reference {input.reference}"
        " --output {output}"
        # " --mask {input.mask}"

rule concentration_estimate:
    input:
        T1files=[
            Path(f"DATA/{subject}/CONCENTRATIONS/{file.with_suffix('.mgz').name}")
            for subject, files in SUBJECT_FILES.items()
            for file in files["T1"]
        ],


rule recon_all_workflow:
    input:
        t1=lambda x: sorted(Path(f"DATA/{x.subjectid}/RESAMPLED").iterdir())[0],
        t2="DATA/{subjectid}/RESAMPLED/T2.mgz",
        script="mri_recon.py",
    output:
        protected(directory("DATA/freesurfer/{subjectid}/")),
    threads: config["recon_threads"]
    resources:
        time="96:00:00",
    params:
        subjectsdir=config["FS_subjectsdir"],
    shadow:
        config["shadow"]
    shell:
        "python {input.script}"
        " --subjectid {wildcards.subjectid}"
        " --t1file {input.t1}"
        " --outputdir {output}"
        " --subjectsdir {params.subjectsdir}"
        # " -parallel"
        # " --t2file {input.t2}"


rule symlink_freesurfer_workflow:
    input:
        "DATA/freesurfer/{subjectid}",
    output:
        directory("DATA/{subjectid}/freesurfer"),
    shell:
        "ln -s '../../{input}' '{output}'"


rule mri_recon_all:
    input:
        expand("DATA/{subjectid}/freesurfer", subjectid=config["subjects"]),


####################################
# Converting data to FEniCS-formats.
####################################

rule mri2fenics:
    input:
        meshfile="DATA/{subjectid}/MODELING/resolution{res}/mesh.hdf",
        timestamps="DATA/{subjectid}/timestamps.txt",
        concentration_files= lambda wc: [
            f"DATA/{wc.subjectid}/CONCENTRATIONS/{file.with_suffix('.mgz').name}"
              for file in T1_FILES[wc.subjectid]
        ],
    output:
        hdf="DATA/{subjectid}/MODELING/resolution{res}/data.hdf",
        visual="DATA/{subjectid}/MODELING/resolution{res}/visual/data.xdmf",
    shell:
        "python gmri2fem/mriprocessing/mri2fenics.py"
        " --mris {input.concentration_files}"
        " --mesh_hdf {input.meshfile}"
        " --output_hdf {output.hdf}"
        " --timestamps {input.timestamps}"

rule mri2fenics_all:
    input:
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/data.hdf",
            subjectid=config["subjects"],
            res=config["resolution"],
        ),


rule surface_stl_conversion:
    input:
        "DATA/{subjectid}/freesurfer/surf/{filename}",
    output:
        "DATA/{subjectid}/MODELING/surfaces/{filename}.stl",
    shell:
        "mris_convert --to-scanner {input} {output}"


rule ventricle_extraction:
    input:
        "DATA/{subjectid}/freesurfer/mri/wmparc.mgz",
    output:
        fs="DATA/{subjectid}/MODELING/surfaces/ventricles",
        stl="DATA/{subjectid}/MODELING/surfaces/ventricles.stl",
    shell:
        "bash scripts/extract-ventricles.sh {input} {output.fs}"
        " && mris_convert --to-scanner {output.fs} {output.stl}"


SURFACES = ["lh.pial", "rh.pial", "lh.white", "rh.white", "ventricles"]


rule create_surfaces:
    input:
        expand(
            "DATA/{subjectid}/MODELING/surfaces/{filename}.stl",
            filename=SURFACES,
            subjectid=config["subjects"],
        ),


rule create_mesh:
    input:
        expand(
            "DATA/{{subjectid}}/MODELING/surfaces/{filename}.stl",
            filename=SURFACES,
        ),
    output:
        "DATA/{subjectid}/MODELING/resolution{res}/mesh.hdf",
    resources:
        time="00:30:00",
    shadow:
        config["shadow"]
    shell:
        "python gmri2fem/mriprocessing/mesh_generation.py"
        " --surfaces {input}"
        " --output {output}"
        " --resolution {wildcards.res}"


rule create_subject_meshes:
    input:
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/mesh.hdf",
            subjectid=config["subjects"],
            res=config["resolution"],
        ),



###########################
# Simulation postprocessing
###########################
rule fenics2mri_workflow:
    input:
        referenceimage="DATA/{subjectid}/freesurfer/mri/T1.mgz",
        simulationfile="DATA/{subjectid}/MODELING/resolution{res}/{funcname}.hdf",
        timestampfile="DATA/{subjectid}/timestamps.txt",
    output:
        "DATA/{subjectid}/MODELING/resolution{res}/MRIs/{funcname}_{idx}.nii.gz",
    shell:
        "python gmri2fem/mriprocessing/fenics2mri.py"
        " --simulationfile {input.simulationfile}"
        " --output {output}"
        " --referenceimage {input.referenceimage}"
        " --timestamps {input.timestampfile}"
        " --timeidx {wildcards.idx}"
        " --functionname 'total_concentration'"


rule fenics2mri:
    input:
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/MRIs/{funcname}_{idx}.nii.gz",
            subjectid=config["subjects"],
            res=config["resolution"],
            funcname="multidiffusion_total_data",
            idx=range(5),
        ),


rule solute_quantification_workflow:
    input:
        "DATA/{subjectid}/MODELING/resolution{res}/{funcname}.hdf",
    output:
        "DATA/{subjectid}/MODELING/resolution{res}/quantities_{funcname}.csv",
    threads: config["sim_threads"]
    shadow:
        config["shadow"]
    shell:
        "OMP_NUM_THREADS=1 mpirun -n {threads}"
        " python gmri2fem/analysis/solute_quantification.py"
        " --input {input} --funcname total_concentration --output {output}"


rule solute_quantification:
    input:
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/quantities_{funcname}.csv",
            funcname=["data", "diffusion", "multidiffusion_total"],
            subjectid=config["subjects"],
            res=config["resolution"],
        ),


rule solute_quantification_data:
    input:
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/quantities_{funcname}_data.csv",
            funcname=["data", "diffusion", "multidiffusion_total", "diffusion_singlecomp"],
            subjectid=config["subjects"],
            res=config["resolution"],
        ),


rule hdf2xdmf:
    input:
        "DATA/{subjectid}/MODELING/resolution{res}/{dataset}.hdf",
    output:
        directory("DATA/{subjectid}/MODELING/resolution{res}/{dataset}_visual/"),
    params:
        compartments="ecs pvs",
        funcname="multidiffusion",
    shell:
        "python3 scripts/hdf2xdmf.py"
        " --input {input}"
        " --outputdir {output}"
        " --funcname {params.funcname}"
        " --subnames {params.compartments}"


rule visualize_differences:
    input:
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/visual/diffusion-data.xdmf",
            subjectid=config["subjects"],
            res=config["resolution"],
        ),
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/visual/multidiffusion_total-data.xdmf",
            subjectid=config["subjects"],
            res=config["resolution"],
        ),
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/visual/multidiffusion_total-diffusion.xdmf",
            subjectid=config["subjects"],
            res=config["resolution"],
        ),


rule visualize_compartment_differences:
    input:
        "DATA/{subjectid}/MODELING/resolution{res}/multidiffusion.hdf",
    output:
        "DATA/{subjectid}/MODELING/resolution{res}/visual/ECS-PVS.xdmf",
    threads: config["sim_threads"]
    shell:
        "python3 gmri2fem/analysis/visualize_differences.py "
        "--input {input} "
        "--funcname 'fluid_concentration'"
        "--fname1 ECS "
        "--fname2 PVS "
        "--output {output}"


# Interpolation error investigation
rule fenics2mri_transpose:
    input:
        expand(
            "DATA/{subjectid}/MODELING/resolution{res}/MRIs/data_{idx}.nii.gz",
            subjectid=config["subjects"],
            res=config["resolution"],
            idx=range(5)
        )

