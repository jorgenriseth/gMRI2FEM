# gMRI2FEM
    
This repository contains code for processing glymphatic MRI-images, i.e. contrast-enhanced brain images, with special focus on concentration-estimation and conversion into dolfin for mathematical modelling of brain tracer transport using the finite element method. 

While this repository may be installed as a python-package and used freely as such, it also contains a `snakemake` processing pipeline built around estimating concentrations from T1-weighted images or T1-maps, and mapping them to `dolfin`-functions. This, of course, relies on having the available data and organizing it according to the data section below.
