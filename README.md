Requirements:
- dcm2niix



## Minor Update

## MRI Preprocessing

All steps are collected into one single script: `mri_preprocessing.py`. Run with
```bash
python mri_preprocessing.py PAT_XXX --skipdicom
```

1. Extract MR-images from DICOM. This is done using the script `multiframe_dicom.py` (assuming the DICOM images are in "multiframe" or "enhanced" format. Otherwise, see the outdated `sort_mri_old.py` for traditional format). The images are extracted in `.nii`-format.

### Concentration Estimation
After the preprocessing steps, ending with registration has been performed, we want to reconstruct concentrations. For this we need to create a refroi-mri using either freeview (according to Bastian's chapter in soon to be available brain-book), or using 3DSlicer. The latter has the advantage of having a flood-fill algorithm making it easy to create a "robust" refroi. The refroi should be stored as `data/PAT_XXX/refroi.nii` (TODO: Make it possible to have both `nii` and `mgz`). This enables normalization:
```bash
python gmri2fem/mriprocessing/normalize_images.py PAT_XXX
```
followed by concentration reconstructions:
```bash
python estimatec.py [collection of arguments to be better determined]
```

## Run Freesurfer recon-all


## Convert to FEniCS-compatible data
Finally, the concentrations should be converted fenics-functions using the `mri2fenics.py`-script.



## Running on a cluster

**Work in progress**

The following instructions are meant to provide an easy guide to run the MRI-processing pipeline using one of the Sigma2-administered clusters such as Saga or Colossus. While some parts only require a FreeSurfer-installation, the later parts of the processing pipeline require various python-packages. We therefor assume installation of these as well.


To activate conda-environment on saga:
```bash
module load Mamba/4.14.0-0
source ${EBROOTMAMBA}/bin/activate
conda activate ../conda/packaga-cache/gmri2fem
echo $EBROOTMAMBA
```

To execute snakemake workflows, we use the `sagamake.sh` script, and can be started as a slurm job in the following way:
```bash
sbatch sagamake.sh path/to/output
```

## Running pipeline using snakemake
### 1. Resample/conform
Assuming the MRI-data has ben extracted from DICOM to nii or mgz, located in 'DATA/{subject_id}/'. Update `snakeconfig.yml`, especially the list of subjects to all that you want to process.
```bash
snakemake mri_convert --cores [ncores]
```
This should run mri_conform on both T1 and T2-data.

### 2. Register
Secondly, perform the registration
```bash
snakemake register_all --cores [ncores]
```

### 3. Create a reference ROI for normalization.
Using freeview, mark a region within the fatty regions behind the eye, to which the images should be normalized. Should be stored at top-level of subject-data, i.e. `DATA/{subject}/refroi.mgz`.

### 4. Normalize
```bash
snakemake mri_normalize --cores [ncores]
```
Inputs registered T1-weighted images `DATA/{subject}/REGISTERED/xxxx.mgz` and the `refroi.mgz`, and outputs the images `DATA/{subject}/NORMALIZED/xxxx.mgz` that are normalized with respect to the median intensity of each image w.r.t. the reference ROI.

### 5. Recon-all
```bash
snakemake mri_recon_all --cores [ncores]
```
Takes as input the first of the T1-weighted images, and potentially the T2-weighted image, and runs recon-all.

The directory will be output to `DATA/freesurfer/{subject}`, and then a symlink will be generated to the directory within each of the subject-folders, i.e. at `DATA/{subject}/freesurfer/ -> DATA/freesurfer/subject`. The reason for this, is that freesurfer assumes a very specific structure of it's data, so this will be useful for e.g. DTI-processing.


### 6. Concentration estimation
```bash
snakemake concentration_estimate --cores [ncores]
```