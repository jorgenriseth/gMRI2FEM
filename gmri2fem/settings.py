from pathlib import Path

from pydantic import BaseModel
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    rawdata: Path = Path(__file__).parent.parent / "GRIP"
    datapath: Path = Path(__file__).parent.parent / "DATA"


class PatientDataSettings(BaseModel):
    dicompath: Path
    datapath: Path
    patient_root: Path
    t1raw: Path
    looklocker: Path
    t2: Path
    resampled: Path
    registered: Path
    normalized: Path
    lta: Path
    dti: Path
    concentrations: Path
    modeling: Path
    freesurfer: Path


def patient_data_default_settings() -> dict[str, str]:
    return dict(
        t1raw="MRIs/T1",
        looklocker="MRIs/LOOKLOCKER",
        t2="MRIs/T2",
        resampled="RESAMPLED",
        registered="REGISTERED",
        normalized="NORMALIZED",
        lta="LTA",
        dti="DTI",
        concentrations="CONCENTRATIONS",
        modeling="MODELING/",
        freesurfer="freesurfer",
    )


def patient_data_settings(patientid: str) -> PatientDataSettings:
    default = patient_data_default_settings()
    patient_path = Settings().datapath / patientid
    return PatientDataSettings(
        dicompath=Settings().rawdata / patientid,
        datapath=Settings().datapath,
        patient_root=patient_path,
        **{data: patient_path / datadir for data, datadir in default.items()}
    )


class DICOMSettings(BaseModel):
    paths: PatientDataSettings
    enchanced: bool = True
    patterns: dict[str, str] = {
        "T1": r"(DICOM_\d+_\d+[_ ])(.*(PDT1_3D).*)",
        "T2": r"(DICOM_\d+_\d+[_ ])(.*(TE565).*)",
        "LookLocker": r"(DICOM_\d+_3[_ ])(.*(LookLocker|2beatpause).*)",
    }
