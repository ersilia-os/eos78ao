import csv
import sys
from rdkit import Chem
from mordred import Calculator, descriptors
from timeout_decorator import timeout
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r") as f:
    reader = csv.reader(f)
    next(reader)
    smiles = []
    for r in reader:
        smiles += [r[0]]

TIMEOUT_SEC=3

calc = Calculator(descriptors, ignore_3D=True)

@timeout(TIMEOUT_SEC)
def one_molecule(mol):
    return calc(mol)


columns = list(calc.pandas([Chem.MolFromSmiles("CCCC")]).columns)

R = []
for smi in smiles:
    try:
        r = one_molecule(Chem.MolFromSmiles(smi))
    except:
        r = [None for _ in range(len(columns))]
    R += [r]

# We read as float dtype because mordred gives us errors for many of the MAX... type descriptors
# Refer to https://github.com/mordred-descriptor/mordred/issues/68
df = pd.DataFrame(R, columns=columns, dtype=float)
# Then we fill those empty columns with 0.0
df.fillna(0.0, inplace=True)
df.to_csv(outfile, index=False)
