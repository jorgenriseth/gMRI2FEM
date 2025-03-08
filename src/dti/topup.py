import tempfile
import re
from pathlib import Path
from typing import Optional

import subprocess

from dti.utils import mri_number_of_frames, readout_time, with_suffix


def topup(dti, topup_b0, outputdir: Path, tmppath: Optional[Path] = None):
    if tmppath is None:
        tmpdir = tempfile.TemporaryDirectory()
        tmppath = Path(tmpdir.name)

    dti_b0 = tmppath / dti.name.replace("DTI.nii.gz", "b0.nii.gz")
    cmd_extract_b0 = f"fslroi {dti} {dti_b0} 150 10"
    subprocess.run(cmd_extract_b0, shell=True, check=True)

    str_base = common_prefix(Path(dti).name, Path(topup_b0).name)
    merged_b0 = tmppath / f"{str_base}_topup_b0_stack.nii.gz"
    cmd = f"fslmerge -t {merged_b0} {dti_b0} {topup_b0}"
    subprocess.run(cmd, shell=True, check=True)

    dti_json = with_suffix(dti, ".json")
    acq_params = outputdir / f"{str_base}_topup_acq_params.txt"
    nframes = mri_number_of_frames(dti_b0)
    main_readout_time = readout_time(dti_json)
    topup_readout_time = readout_time(with_suffix(topup_b0, ".json"))
    acq_params.parent.mkdir(exist_ok=True, parents=True)
    with open(acq_params, "w") as f:
        for i in range(nframes):
            f.write(f"0 -1 0 {main_readout_time}\n")
        f.write(f"0 1 0 {topup_readout_time}\n")

    output = outputdir / f"{str_base}_topup_output"
    topup_cmd = (
        f"topup"
        + f" --imain={merged_b0}"
        + f" --datain={acq_params}"
        + f" --config=b02b0.cnf"
        + f" --out={output}"
        + f" --iout={with_suffix(output, '_b0')}"
        + f" --fout={with_suffix(output, '_field.nii.gz')}"
    )
    subprocess.run(topup_cmd, shell=True, check=True)


def common_prefix(str1, str2):
    for idx, (c1, c2) in enumerate(zip(str1, str2)):
        if c1 != c2:
            parts = str1[:idx].split("_")
            if re.match(r"\w+-\w*", parts[-1]):
                return "_".join(parts[:-1])
            return str1[:idx]
    return ""
