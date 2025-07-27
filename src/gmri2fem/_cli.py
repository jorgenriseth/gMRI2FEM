# ruff: disable=F401

import click

from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "mixed-t1map": "gmri2fem.mixed_t1map.mixed_t1map",
        "looklocker-t1map": "gmri2fem.looklocker_t1map.looklocker_t1map",
        "hybrid-t1map": "gmri2fem.hybrid_t1map.hybrid_t1map",
        "looklocker-t1-postprocessing": "gmri2fem.t1maps.looklocker_t1_postprocessing",
        "t1-to-r1": "gmri2fem.t1maps.T1_to_R1",
        "t1w-sigdiff": "gmri2fem.t1_weighted.T1w_sigdiff",
        "t1w-normalize": "gmri2fem.t1_weighted.T1w_normalize",
        "concentration": "gmri2fem.concentration.concentration",
        "reslice4d": "gmri2fem.reslice_4d.reslice4d",
        "stats": "gmri2fem.regionwise_statistics.compute_mri_stats",
    },
)
def mri():
    pass


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "refine": "gmri2fem.segmentation_refinement.refine",
        "mask-intracranial": "gmri2fem.masking.mask_intracranial",
        "mask-csf": "gmri2fem.masking.mask_csf",
        "orbital-refroi": "gmri2fem.orbital_refroi.orbital_refroi",
        "extended-fs": "gmri2fem.csf_segmentation.create_extended_segmentation_cli",
    },
)
def seg():
    pass
