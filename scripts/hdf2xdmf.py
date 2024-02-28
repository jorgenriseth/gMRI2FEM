import argparse
from pathlib import Path

from pantarei.fenicsstorage import FenicsStorage

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True)
parser.add_argument("--outputdir", type=str, required=True)
parser.add_argument("--funcname", type=str, required=True)
parser.add_argument("--subnames", type=str, nargs="+", required=True)
args = parser.parse_args()

file = FenicsStorage(Path(args.input), "r")
file.to_xdmf(args.funcname, args.subnames, outputpattern=Path(args.outputdir))
file.close()
