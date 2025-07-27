CORPUS_CALLOSUM = [
    251,
    252,
    253,
    254,
    255,
]
CEREBRAL_WM_RANGES = [
    *[2, 41],  # aseg left/right cerebral white labels
    *list(range(3000, 3036)),  # wmparc-left-labels
    *list(range(4000, 4036)),  # wmparc-right-labels
    *[5001, 5002],
    *[28, 60],  # VentralDC included in white matter subdomain
    *CORPUS_CALLOSUM,  #  Corpus callosum
    *[31, 63],  # Choroid plexus.
]
CEREBRAL_CGM_RANGES = [
    *[3, 42],  # aseg left/right cortcal gm
    *list(range(1000, 1036)),  # aparc left labels
    *list(range(2000, 2036)),  # aparc right labels
]
SUBCORTICAL_GM_RANGES = [
    *(10, 49),  # Thalamus,
    *(11, 50),  # Caudate,
    *(12, 51),  # Putamen,
    *(13, 52),  # pallidum
    *(17, 53),  # hippocampus
    *(18, 54),  # amygdala
    *(26, 58),  # accumbens
]
VENTRICLES = [
    4,  # Left lateral ventricle
    5,  # Left inferior lateral ventricle
    14,  # Third ventricle
    15,  # Fourth ventricle
    43,  # Right lateral ventricle
    44,  # Right inferior lateral ventricle
]
FREESURFER_CSF = [
    *VENTRICLES,
    24,  # Generic CSF
]
CORTICAL_CSF = [
    *(x + 15000 for x in CEREBRAL_CGM_RANGES),
]
SEGMENTATION_GROUPS = {
    "cerebral-wm": CEREBRAL_WM_RANGES,
    "cerebral-cortex": CEREBRAL_CGM_RANGES,
    "cerebellar-wm": [7, 46],  # left, right
    "cerebellar-cortex": [8, 47],  # left, right
    "csf-freesurfer": FREESURFER_CSF,
    "cortical-csf": CORTICAL_CSF,
    "corpus-callosum": CORPUS_CALLOSUM,
}


def default_segmentation_groups():
    groups = {**SEGMENTATION_GROUPS}
    return groups


# Not certain if useful, kept until first commit
def invert_seg_groups(seg_groups):
    d = {}
    for key, vals in seg_groups.items():
        for val in vals:
            if val not in d:
                d[val] = [key]
            else:
                d[val].append(key)
    return d


if __name__ == "__main__":
    import json

    print(json.dumps(default_segmentation_groups(), indent=4))
