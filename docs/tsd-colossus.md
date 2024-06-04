
# Setting up gMRI2FEM for TSD
## Installing snakemake through conda
1. Install `conda` in the user home-directory on TSD by importing Miniforge-installer into TSD.
2. On personal computer, create a `conda`-environment with python<3.12 and snakemake7.x.x.
3. Package environment into a `tar.gz` using `conda-pack`:
```bash
mamba pack -n env_name -o your_env_packed.tar.gz
```
4. Import tar-file, and unpack it into a directory
```bash
cd [.../]miniforge3/envs/
mkdir env_name
tar xvfz yourenv.tar.gz -C env_name
```
The environment should now be possible to activate through `conda`/`mamba`.

## Packaging and importing files
### Dataset folder
Since there might be variations in the DICOM source data, the `mri_dataset` folder is expected to contain its own `Snakefile`, and code for extracting DICOM data and immediate derivatives. Necessary files to import:
- `code/`
- `Snakefile`
```bash
tar cvfz mri_dataset_code.tar.gz code/ Snakefile
```
Together with a `participants-private.json`-file generated from the overview of contrast injection times, the `Snakefile` and the scripts should hopefully suffice to generate necessary metadata, timetables etc.

### GMRI2FEM
1. Find exact list of files that is required to run
```bash
tar cvfz gmri2fem.tar.gz \
  gmri2fem \ 
  scripts \ 
  snakeprofiles \ 
  singularity/{gmri2fem.sif,license.txt} \ 
  environment.yml \ 
  pyproject.toml \ 
  sagamake.sh \ 
  Snakefile \ 
  snakeconfig.yaml \ 
  workflows_additional
```

## Running workflows on TSD 
```bash
ssh pXXXX-submit.tsd.usit.no
```
```bash
source .bashrc
conda activate [snakemake_env]
```
It is now possible to run `snakemake` from the login-node for small jobs or dry-runs. For larger jobs, we recommend using the HPC-cluster `colossus`. To run slurm-jobs using snakemake, adjust the needed resources in `snakeprofiles/colossus/config.yaml`, and execute the desired workflows by replacing `snakemake` by `sbatch sbatches/colossusmake.sh`. Additional `snakemake` arguments may be provided as usual.

```bash
sbatch sbatches/colossusmake.sh [--config subjects=["sub-XX", "sub-XY"]]
```
