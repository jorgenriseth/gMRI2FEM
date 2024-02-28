import argparse
from pathlib import Path

import numpy as np

SUBDOMAIN_LABELS = {
    "gray-matter": 1,
    "white-matter": 2,
}

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, required=True)


from dataclasses import dataclass, field


@dataclass
class DummyArgs:
    input: Path
    output: Path


args = DummyArgs(
    input=Path("DATA/CTRL_001/MODELING/resolution16/data.hdf"),
    output=Path("DATA/CTRL_001/MODELING/resolution32/fem_quantities.json"),
)

import dolfin as df
import pantarei as pr

# %%
SUBDOMAIN_LABELS = {
    "gray-matter": 1,
    "white-matter": 2,
}


def fem_solute_quantities(
    input: Path, subdomain_labels: dict[str, int]
) -> tuple[list[float], dict[str, dict[str, int | list[int] | float | list[float]]]]:
    with df.HDF5File(df.MPI.comm_world, str(input), "r") as hdf:
        domain = pr.read_domain(hdf)
        dx = df.Measure("dx", domain, subdomain_data=domain.subdomains)
        timevec = pr.read_timevector(hdf, "total_concentration")
        u = pr.read_function(hdf, "total_concentration", domain)
        region_data = {
            region: {
                "idx": region_idx,
                "volume": df.assemble(1 * dx(region_idx)),
                "mass": np.nan * np.ones_like(timevec),
                "mean": np.nan * np.ones_like(timevec),
            }
            for region, region_idx in subdomain_labels.items()
        }
        for idx, ti in enumerate(timevec):
            pr.read_checkpoint(hdf, u, "total_concentration", idx)
            for region, data in region_data.items():
                data["mass"][idx] = df.assemble(u * dx(data["idx"]))
                data["mean"][idx] = data["mass"][idx] / data["volume"]

    region_data["all"] = {
        "idx": list(data["idx"] for data in region_data.values()),
        "volume": sum((data["volume"] for data in region_data.values())),
        "mass": sum((data["mass"] for data in region_data.values())),
        "mean": sum((data["mean"] for data in region_data.values())),
    }
    return timevec, region_data


# %%
import matplotlib.pyplot as plt

time, region_data = fem_solute_quantities(args.input, SUBDOMAIN_LABELS)
plt.figure()
for region, data in region_data.items():
    plt.plot(time, data["mean"], label=region)
plt.legend()
plt.show()

import json

with open(args.output, "w") as f:
    region_data["time"] = list(time)
    json.dump(region_data, f, indent=4)
