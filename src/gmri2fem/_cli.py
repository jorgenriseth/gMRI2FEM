# ruff: disable=F401
import click
from pathlib import Path
from typing import Optional
from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "mixed-t1map": "gmri2fem.mixed_t1map.mixed_t1map",
        "looklocker-t1map": "gmri2fem.looklocker_t1map.looklocker_t1map",
        "hybrid-t1map": "gmri2fem.hybrid_t1map.hybrid_t1map",
        "looklocker_t1_postprocessing": "gmri2fem.t1maps.looklocker_t1_postprocessing",
        "t1-to-r1": "gmri2fem.t1maps.T1_to_R1",
        "t1w-sigdiff": "gmri2fem.t1_weighted.T1w_sigdiff",
        "t1w-normalize": "gmri2fem.t1_weighted.T1w_normalize",
        "concentration": "gmri2fem.concentration.concentration",
        "reslice4d": "gmri2fem.reslice_4d.reslice4d",
    },
)
def mri():
    from gmri2fem.looklocker_t1map import looklocker_t1map
    from gmri2fem.mixed_t1map import mixed_t1map
    from gmri2fem.hybrid_t1map import hybrid_t1map
    from gmri2fem.t1maps import T1_to_R1, looklocker_t1_postprocessing
    from gmri2fem.t1_weighted import T1w_sigdiff, T1w_normalize
    from gmri2fem.concentration import concentration
    from gmri2fem.reslice_4d import reslice4d

    mri.add_command(mixed_t1map)
    mri.add_command(looklocker_t1map)
    mri.add_command(looklocker_t1_postprocessing)
    mri.add_command(hybrid_t1map)
    mri.add_command(T1_to_R1)

    mri.add_command(T1w_sigdiff)
    mri.add_command(T1w_normalize)

    mri.add_command(concentration)
    mri.add_command(reslice4d)


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "refine": "gmri2fem.segmentation_refinement.refine",
        "mask-intracranial": "gmri2fem.masking.mask_intracranial",
        "mask-csf": "gmri2fem.masking.mask_csf",
        "orbital-refroi": "gmri2fem.orbital_refroi.orbital_refroi",
    },
)
def seg():
    from gmri2fem.masking import mask_intracranial, mask_csf
    from gmri2fem.segmentation_refinement import refine
    from gmri2fem.orbital_refroi import orbital_refroi

    seg.add_command(refine)
    seg.add_command(mask_intracranial)
    seg.add_command(mask_csf)
    seg.add_command(orbital_refroi)


@click.command()
@click.argument("dcmpath", type=Path, required=True)
@click.argument("outpath", type=Path, required=True)
@click.option("--subvolumes", type=str, multiple=True, default=None)
def dcm2nii_mixed(
    dcmpath: Path,
    outpath: Path,
    subvolumes: Optional[list[str]] = None,
):
    import json
    import subprocess
    import nibabel
    from gmri2fem.mixed_dicom import VOLUME_LABELS, dcm2nii_mixed

    subvolumes = subvolumes or VOLUME_LABELS
    assert all(
        [volname in VOLUME_LABELS for volname in subvolumes]
    ), f"Invalid subvolume name in {subvolumes}, not in {VOLUME_LABELS}"
    outdir, form = outpath.parent, outpath.stem
    outdir.mkdir(exist_ok=True, parents=True)

    vols = dcm2nii_mixed(dcmpath, subvolumes)
    meta = {}
    for vol, volname in zip(vols, subvolumes):
        output = outpath.with_name(outpath.stem + "_" + volname + ".nii.gz")

        nii = vol["nifti"]
        descrip = vol["descrip"]
        nibabel.nifti1.save(nii, output)
        if volname == "SE-modulus":
            meta["TR_SE"] = descrip["TR"]
            meta["TE"] = descrip["TE"]
        elif volname == "IR-corrected-real":
            meta["TR_IR"] = descrip["TR"]
            meta["TI"] = descrip["TI"]

    with open(outpath.parent / f"{form}_meta.json", "w") as f:
        json.dump(meta, f)

    try:
        cmd = f"dcm2niix -w 0 --terse -b o -f '{form}' -o '{outdir}' '{dcmpath}' >> /tmp/dcm2niix.txt "
        subprocess.run(cmd, shell=True).check_returncode()
    except (ValueError, subprocess.CalledProcessError) as e:
        print(str(e))
        pass

    pass
