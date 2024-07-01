# Running containers on clusters
Several of the workflows requires large amounts of memory, and runs for a long time. It is therefore preferable to defer the work to a cluster, to keep your personal laptop/workstation available for use.

## Prepare the files to be sent 
To run the workflow, we are going to need a few files, that aren't necessarily available on the cluster. For simplicity, define the following env-variables:
```bash
export GMRI2FEM_REMOTE_PATH="[host]:[/path/to/project/in/remote]"
```
Filelist:
```bash
- gmri2fem/
- scripts/
- snakeprofiles/
- singularity/gmri2fem.sif
- singularity/license.txt
- Snakefile
- snakeconfig.yaml
- envionment.yml
- pyproject.toml
- sagamake.sh
- workflows_additional
```
Send all necessary source files:
```bash
rsync --info=progress2 -za --relative --mkpath gmri2fem scripts snakeprofiles singularity/{gmri2fem.sif,license.txt} environment.yml pyproject.toml sbatches/sagamake.sh Snakefile snakeconfig.yaml workflows_additional $GMRI2FEM_REMOTE_PATH
```
Decide which parts of the data are necessary to send, and send those along as well, using same command as above. Some that are easily forgotten 
```
  tar cvfz data.tar.gz data/mri_dataset/{sub-*,derivatives,participants.json,participants.tsv,timetable.csv,Snakefile}
  rsync  -v --info=progress2 data.tar.gz $GMRI2FEM_REMOTE_PATH
```
Log into remote over ssh, and enter the project directory. Since it's preferable to have a distinct job for executing `snakemake`, we run the commands
```bash
sbatch sagamake.sh [file or rule]
```
where `[file or rule]` refers to a typical target for a snakemake workflow. This will schedule a job which executes `snakemake` with the `saga` snakeprofile, which specifies several command line arguments for starting the jobs with singularity and slurm as well as the available resources.

# Easily move contents of a folder to HPC cluster
Assume that you have a directory `projectName` with all the files etc. needed to run a set of scripts , that you want to get running on an HPC cluster with `snakemake`.

## Quick setup of `snakemake` for clusters.
While we will be aiming to mainly use containers for the `conda`-environments used as part of the workflow, we still require that `snakemake` is available from the login-node. If this is not already avaiable through e.g. `module avail snakemake`, it can be installed through `conda`, which is often available, or can be installed.
While it is often discouraged to install conda into the users home-directory, since the home-directory storage is typically very small and is not as durable as the project-folder, we only require `snakemake` and will therefore ignore the warning, as the systemwide installation is often cumbersome and slow.
```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p $HOME/conda
source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
```



## Executing snakemake workflows with symlinks
When snakemake executes workflows using singularity, it passes the argument "--home [path/to/project]" to `singularity exec`. This means that files whose absolute path is not situated within the proejct root directory that are used by workflows needs to be explicitly passed to 
```bash
snakemake [workflow or filepaths]--cores N --use-singularity --singularity-args "\-c \-B license.txt:/license.txt \-B '$(pwd)' [\-B /proper/path/to/data]"
```
NB: The proper path to data is included in case the data-directory is symlinked. However, if you have a symlink-chain, then all intermediate paths needs to be included. Therefore it is recommended to symlink directly to the proper target. This can be found by `readlink -f [symlink]`
