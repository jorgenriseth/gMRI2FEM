from pathlib import Path
from typing import Optional

import click
import dolfin as df
import fenics_io as fio


def collect_mesh_data(
    domain_data: Path,
    concentration_data: Path,
    parcellation_data: Path,
    output: Path,
    dti_data: Optional[Path] = None,
):
    comm = df.MPI.comm_world
    domain = fio.read_domain(domain_data)
    with df.HDF5File(comm, str(output), "w") as hdf:
        fio.write_domain(hdf, domain)

        with df.HDF5File(comm, str(concentration_data), "r") as concentration_hdf:
            for funcname in ["concentration", "boundary_concentration"]:
                func = fio.read_function(concentration_hdf, funcname, domain)
                time = fio.read_timevector(concentration_hdf, funcname)
                fio.write_function(hdf, func, funcname)
                for idx, ti in enumerate(time[1:], start=1):
                    fio.read_checkpoint(concentration_hdf, func, funcname, idx)
                    fio.write_checkpoint(hdf, func, funcname, ti)

        if dti_data is not None:
            with df.HDF5File(comm, str(dti_data), "r") as dti_hdf:
                for funcname in ["DTI", "FA", "MD"]:
                    func = fio.read_function(dti_hdf, funcname, domain)
                    fio.write_function(hdf, func, funcname)

        with df.HDF5File(comm, str(parcellation_data), "r") as parcellation_hdf:
            subdomains = df.MeshFunction("size_t", domain, domain.topology().dim())
            parcellation_hdf.read(subdomains, "parcellations")
        hdf.write(subdomains, "parcellations")


@click.command()
@click.option("--domain", "domain_data", type=Path, required=True)
@click.option("--concentration_data", type=Path, required=True)
@click.option("--parcellation_data", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--dti_data", type=Path)
def collect(*args, **kwargs):
    collect_mesh_data(*args, **kwargs)
