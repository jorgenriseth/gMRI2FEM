import csv
import json
from pathlib import Path
from typing import Any


def participants_header(filepath: Path) -> list[str]:
    with open(filepath, "r") as f:
        d = json.load(f)
        header = ["subject-id", *list(d.keys())]
    return header


def write_records_to_tsv(
    records: list[dict[str, Any]], header: list[str], filepath: Path
):
    with open(filepath, "w") as f:
        w = csv.DictWriter(f, header, delimiter="\t")
        w.writeheader()
        for record in records:
            w.writerow(record)
    return filepath


def generate_metadata(
    participants_private: Path,
    participants_json: Path,
    participants_tsv: Path,
):
    with open(participants_private, "r") as f:
        private_data = json.load(f)
    header = participants_header(participants_json)
    records = [{key: subject[key] for key in header} for subject in private_data]
    write_records_to_tsv(records, header, participants_tsv)
    return private_data


def read_records_from_tsv(filepath: Path) -> list[dict[str, str]]:
    with open(filepath, "r") as f:
        header = next(csv.reader(f, delimiter="\t"))
        r = csv.DictReader(f, fieldnames=header, delimiter="\t")
        records = [record for record in r]
    return records


if __name__ == "__main__":
    generate_metadata(
        participants_private=Path("participants-private.json"),
        participants_json=Path("participants.json"),
        participants_tsv=Path("participants.tsv"),
    )
