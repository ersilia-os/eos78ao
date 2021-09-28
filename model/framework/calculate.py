import csv
import sys
from rdkit import Chem
from mordred import Calculator, descriptors

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r") as f:
    reader = csv.reader(f)
    next(reader)
    smiles = []
    for r in reader:
        smiles += [r[0]]

calc = Calculator(descriptors, ignore_3D=True)

df = calc.pandas([Chem.MolFromSmiles(smi) for smi in smiles])

df.to_csv(outfile, index=False)
