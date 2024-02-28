#!/usr/bin/env bash
# System setup
# set -o errexit  # Exit the script on any error
# set -o nounset  # Treat any unset variables as an error
set -euf -o pipefail # settings to catch errors in bash scripts

# Define paths
PROJECTSDIR="/cluster/projects/($PROJECTID)/($USER)"
PROJECTNAME="gmri2fem"

module load Mamba/4.14.0-0
#module load Minicoda3/22.11.1-1
export CONDA_PKGS_DIRS=${PROJECTSDIR}/conda/packaga-cache/gmri2fem
export CONDA_ENVS_DIRS=${PROJECTSDIR}/conda

echo ${CONDA_ENVS_DIRS}/${PROJECTNAME}
mamba env create --prefix ${CONDA_PKGS_DIRS}/${PROJECTNAME} --file environment.yml # <-- used for creating the project at the beginning...
#mamba env update --prefix ${CONDA_STORAGE}/${PROJECTNAME} --file environment.yml  # <-- ... but only needs updating from env-file later.
