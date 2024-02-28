import os
import re
from pathlib import Path

import pandas as pd
import dataclasses


@dataclasses.dataclass
class LutEntry:
    index: int
    description: str


def find_lut_entry(line: str) -> LutEntry | None:
    m = re.match(r"^(\d+)\s+([_\da-zA-Z-]+)\s.", line)
    if m is None:
        return None
    return LutEntry(index=int(m.group(1)), description=m.group(2))


def to_dict(dataframe: pd.DataFrame) -> dict[int, str]:
    return dataframe.set_index("lut_index").to_dict()["description"]


if __name__ == "__main__":
    freesurfer_lut = Path(os.environ["FREESURFER_HOME"]) / "FreeSurferColorLUT.txt"
    with open(freesurfer_lut, "r") as f:
        lut_entries = [x for x in map(find_lut_entry, f) if x is not None]
        indices = pd.Series(map(lambda x: x.index, lut_entries), dtype=int)
        descriptions = pd.Series(map(lambda x: x.description, lut_entries), dtype=str)

    df = pd.DataFrame({"lut_index": indices, "description": descriptions})
    df.set_index("lut_index")["description"].to_json(
        "data/freesurfer_lut.json", indent=4
    )
    df.to_csv("data/freesurfer_lut.csv", index=False)
