# gMRI2FEM
`gmri2fem` is a python package for processing of human glymphatic MRI-images, i.e. contrast-enhanced brain images, with special focus on concentration-estimation and conversion into `dolfin` for mathematical modelling of brain tracer transport using the finite element method.
The data functionality was developed for the MRI processing workflow for the GRIP-project. 
The package both provide a library of useful functions as well as a CLI for processing different types of images.


## Installation
### Dependencies
`gmri2fem` has some python- and non-python dependencies which are not easily installable through `pip`:
- FreeSurfer
- greedy
- FEniCS
- SVMTK

Instructions for how to install these may be found on their websites. 
It is, however, possible to `gmri2fem` and run large portions of the pipeline whithout these dependencies.
To install the python package gmri2fem, clone this repository and run 
```bash
pip install . 
```
from the root directory. 
This will install the python packages and the CLI. 


## Example usage
For a rough documentation on how to use the various exposed CLI-commands, we recommend looking at the Snakefiles for the Gonzo data pipeline at (https://github.com/jorgenriseth/gonzo).
The CLI has several levels of subcommands.
To run the main entrypoint:
```bash
$ gmri2fem --help
Usage: gmri2fem [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  brainmeshing
  dcm2nii-ll
  dcm2nii-mixed
  dti
  i2m
  mri
  seg
```

#### MRI
```bash
Usage: gmri2fem mri [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  concentration
  hybrid-t1map
  looklocker-t1-postprocessing
  looklocker-t1map
  mixed-t1map
  reslice4d
  t1-to-r1
  t1w-normalize
  t1w-sigdiff
```

Example command:
```bash
Usage: gmri2fem mri looklocker-t1map [OPTIONS]

Options:
  --input PATH       [required]
  --timestamps PATH  [required]
  --output PATH      [required]
  --help             Show this message and exit.
```