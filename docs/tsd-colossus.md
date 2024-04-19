# Setting up

## Installing snakemake through conda
1. Install conda locally by importing Miniforge installer.
2. On personal computer, create a conda environment with python<3.12 and snakemake7.x.x.
3. Package environment into a tar.gz using conda-pack
```bash
mamba pack -n env_name -o your_env_packed.tar.gz
```
4. Import tar-file, and unpack it into a directory
```bash
cd [.../]miniforge3/envs/
mkdir env_name
tar xvfz yourenv.tar.gz -C env_name
```
It should now be available through conda/mamba

## Packaging and importing files
1. Find exact list of files that is required to run

## Running snakemake 
1. Print necessary command with necessary binds to get it running.
