# Running containers on clusters
Several of the workflows requires large amounts of memory, and runs for a long time. It is therefore preferable to defer the work to a cluster such as TSD-colossus, to keep your personal laptop/workstation available for use.

## Prepare the files to be sent 

To run the workflow, we are going to need a few files, that aren't necessarily available on the cluster. For simplicity, define the following env-variables:
```bash
export GMRI2FEM_REMOTE_PATH="[host]:[/path/to/project/in/remote]"

```
1. Singularity-container: from project root folder, run (also needs license)
  ```bash
  rsync --mkpath --info=progress2 singularity/gmri2fem.sif saga:[/path/to/project/in/remote]/singularity/gmri2fem.sif
  rsync --relative singularity/license.txt $GMRI2FEM_REMOTE_PATH

  ```
2. Required input-data. This two-step procedure might not be the most efficient when `rsync` is available, but is closer to the corresponding procedure on secure clusters without an internet connection.
  ```bash
  tar cvfz data.tar.gz data/mri_dataset/[somesubject] data/mri_processed_data/freesurfer/[somesubject]
  rsync --relative -r data.tar.gz $GMRI2FEM_REMOTE_PATH
  rsync --relative -r data/mri_dataset/derivatives/timetable.csv $GMRI2FEM_REMOTE_PATH

3. Snakefile and neccessary additonals:
  ```bash
  rsync --relative -r Snakefile snakeconfig.yaml snakemake/profiles/[profile] $GMRI2FEM_REMOTE_PATH
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
```
Send all necessary source files:
```bash
rsync --info=progress2 -za --relative --mkpath gmri2fem scripts snakeprofiles singularity/{gmri2fem.sif,license.txt} environment.yml pyproject.toml sagamake.sh Snakefile snakeconfig.yaml $GMRI2FEM_REMOTE_PATH
```
Decide which parts of the data are necessary to send, and send those along as well, using same command as above. Some that are easily forgotten 
```bash
data/mri_dataset/derivatives/timetable.csv
data/mri_dataset/participants.tsv
data/mri_processed_data/freesurfer/[subject] 
```
Log into remote over ssh, and enter the project directory. Since it's preferable to have a distinct job for executing `snakemake`, we run the commands
```bash
sbatch sagamake.sh [file or rule]
```
where `[file or rule]` refers to a typical target for a snakemake workflow. This will schedule a job which executes `snakemake` with the `saga` snakeprofile, which specifies several command line arguments for starting the jobs with singularity and slurm as well as the available resources.


## Singularity/Apptainer
The Sigma2-administered clusters rely on Apptainers (previously Singularity) or podman for running containers. We will be using apptainers. The official documentation can be a handful to get started, but a nice tutorial is available at [https://carpentries-incubator.github.io/singularity-introduction/]

### Convert Dockerfile to Apptainers
As Docker is more broadly used, a lot of the time we will try to build the Apptainers based on Dockerfiles rather than creating the recipe from scratch. 

One solution is to push the desired docker container to Docker Hub, and build it from there using apptainer commands. If you rather want to keep the container local as well, it's possible to save/archive docker-containers to tar-files

Yes, it is possible to build an Apptainer (Singularity) container from a Docker image without pushing the image to Docker Hub or any other remote registry. You can do this by saving the Docker image as a tar archive, transferring the archive to the system where you want to build the Apptainer container, and then using Apptainer to build the container from the tar archive. Here's how you can do it:

1. **Build the Docker Image**: First, build your Docker image from the Dockerfile as you normally would.

   ```bash
   docker build -t your-image-name .
   ```

2. **Save the Docker Image to a Tar Archive**: After building your Docker image, you can save it to a tar file using the `docker save` command.

   ```bash
   docker save your-image-name -o your-image-name.tar.gz | tqdm --bytes --total $(docker image inspect your-image-name --format='{{.Size}}') > your-image-name.tar.gz
   ```
   This command will create a tar archive named `your-image-name.tar` containing your Docker image. The pipe to tqdm is to show a progress bar and can be omitted.

3. **Transfer the Tar Archive**: If the system where you intend to build the Apptainer image is different from the one where you built the Docker image, you'll need to transfer the tar archive to the target system. You can use tools like `scp`, `rsync`, or any other file transfer method.

4. **Build the Apptainer Image from the Tar Archive**: On the target system where Apptainer is installed, use the following command to build the Apptainer image from the tar archive:

   ```bash
   apptainer build your-container-name.sif docker-archive://your-image-name.tar.gz
   ```

   This command will create an Apptainer image named `your-container-name.sif` from the tar archive.

By using this method, you can avoid the need to push your Docker image to a remote registry, which can be particularly useful for development or testing environments, or when working with sensitive or private images that you don't want to upload to a public registry.


# Mamba container requires some stuff to activate upon starting.
```bash
source /opt/conda/etc/profile.d/mamba.sh
source activate [yourenv]
```

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


## Sending files
Use rsync to move the entire directory:
```bash
rsync -az --info=progress2 gmri2fem [host:][host/project/path]
```

## Executing snakemake workflows with symlinks
When snakemake executes workflows using singularity, it passes the argument "--home [path/to/project]" to `singularity exec`. This means that files whose absolute path is not situated within the proejct root directory that are used by workflows needs to be explicitly passed to 
```bash
snakemake [workflow or filepaths]--cores N --use-singularity --singularity-args "-c -B license.txt:/license.txt -B '$(pwd)' [-B /proper/path/to/data]"
```

