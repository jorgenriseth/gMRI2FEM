import numbers
from typing import Any, Optional

import pint

# ureg = pint.UnitRegistry()
ureg = pint.get_application_registry()
dimensionless = pint.Unit("")
mm = ureg.mm
cm = ureg.cm
s = ureg.second
minute = ureg.minute


MULTIDIFFUSION_PARAMETERS = {
    "D_e": 1.3e-4 * mm**2 / s,
    "D_p": 3 * 1.3e-4 * mm**2 / s,
    "phi_e": 0.20 * dimensionless,
    "phi_p": 0.02 * dimensionless,
    "t_ep": 2.9e-2 * 1 / s,
    "t_pb": 0.2e-5 * 1 / s,
    "k_p":  3.7e-4* mm / s,
    "k_e": 1.0e-5 * mm / s,
}

def multidiffusion_parameters(param_units: Optional[dict[str, str]] = None):
    defaults = {**MULTIDIFFUSION_PARAMETERS}
    coefficients = {
        "phi": {
            "ecs": defaults["phi_e"],
            "pvs": defaults["phi_p"]
        },
        "D": {
            "ecs":  defaults["D_e"],
            "pvs": defaults["D_p"],
        },
        "r": {
            "ecs": 0.0 * (1 / s),
            "pvs": defaults["t_pb"],
        },
        "beta": defaults["t_ep"],
        "robin": {
            "ecs": defaults["k_e"],
            "pvs": defaults["k_p"],
        }
    }
    if param_units is not None:
        return make_dimless(convert_to_units(coefficients, param_units))
    return coefficients


def singlecomp_parameters(param_units: Optional[dict[str, str]] = None):
    defaults = {**MULTIDIFFUSION_PARAMETERS}
    coefficients = {
        "D": defaults["D_e"],
        "r": defaults["t_pb"] / defaults["phi_e"],
        "robin": defaults["k_e"] / defaults["phi_e"]
    }
    if param_units is not None:
        return {
            key: val.to(param_units[key]).magnitude
            for key, val in coefficients.items()
        }
    return coefficients


def fasttransfer_parameters(param_units: Optional[dict[str, str]] = None):
    defaults = {**MULTIDIFFUSION_PARAMETERS}
    D_e, D_p = defaults["D_e"], defaults["D_p"]
    phi_e, phi_p = defaults["phi_e"], defaults["phi_p"]
    coefficients = {
        "D": (phi_e * D_e + phi_p * D_p) / (phi_e + phi_p),
        "r": defaults["t_pb"] / (phi_e + phi_p),
        "robin": (defaults["k_p"] + defaults["k_e"]) / (phi_e + phi_p)
    }
    if param_units is not None:
        return {
            key: val.to(param_units[key]).magnitude
            for key, val in coefficients.items()
        }
    return coefficients


def make_dimless(params):
    """Converts all quantities to a dimless number."""
    dimless = {}
    for key, val in params.items():
        if isinstance(val, dict):
            dimless[key] = make_dimless(val)
        else:
            dimless[key] = val.magnitude
    return dimless


def make_quantity_dimless(val: pint.Quantity) -> float:
    assert val.unitless, f"Parameter {val} not dimless"
    return (val.to_base_units()).magnitude


def convert_to_units(params, param_units):
    """Converts all quantities to the units specified by
    param_units."""
    converted = {}
    for key, val in params.items():
        if isinstance(val, dict):
            converted[key] = {}
            for j in val:
                converted[key][j] = val[j].to(param_units[key])
        elif isinstance(val, pint.Quantity):
            converted[key] = val.to(param_units[key])
        else:
            converted[key] = val
    return converted


def isquantity(x):
    return isinstance(x, pint.Quantity) or isinstance(x, numbers.Complex)


def print_quantities(p, offset, depth=0):
    """Pretty printing of dictionaries filled with pint.Quantities"""
    format_size = offset - depth * 2
    for key, value in p.items():
        if isinstance(value, dict):
            print(f"{depth*'  '}{str(key)}")
            print_quantities(value, offset, depth=depth + 1)
        else:
            if isquantity(value):
                print(f"{depth*'  '}{str(key):<{format_size+1}}: {value:.3e}")
            else:
                print(f"{depth*'  '}{str(key):<{format_size+1}}: {value}")


if __name__ == "__main__":
    import pprint
    print("=== Multidiffusion ===")
    print_quantities(multidiffusion_parameters(), 4)

    print()
    print("=== Diffusion - Single compartment===")
    print_quantities(singlecomp_parameters(), 4)

    print()
    print("=== Diffusion - Fast transfer ===")
    print_quantities(fasttransfer_parameters(), offset=4)

