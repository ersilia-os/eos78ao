import os
import csv
import sys
from rdkit import Chem
import pandas as pd
import numpy as np
from mordred import Calculator, descriptors
from timeout_decorator import timeout
import joblib

infile = sys.argv[1]
outfile = sys.argv[2]

ROOT = os.path.abspath(os.path.dirname(__file__))
checkpoints_dir = os.path.join(ROOT, "..","..","checkpoints")

with open(infile, "r") as f:
    reader = csv.reader(f)
    next(reader)
    smiles = []
    for r in reader:
        smiles += [r[0]]

TIMEOUT_SEC=60

calc = Calculator(descriptors, ignore_3D=True)

@timeout(TIMEOUT_SEC)
def one_molecule(mol):
    return calc(mol)

def convert_to_float(df):
    for index, row in df.iterrows():
        for col in df.columns:
            value = row[col]
            if isinstance(value, str):
                try:
                    df.at[index, col] = float(value)
                except ValueError:
                    df.at[index, col] = np.nan
    df = df.map(lambda x: np.nan if pd.isna(x) else x)
    return df

columns = list(calc.pandas([Chem.MolFromSmiles("CCCC")]).columns)

R = []
invalid_idxs = []
for idx, smi in enumerate(smiles):
    try:
        r = one_molecule(Chem.MolFromSmiles(smi))
    except:
        r = [None for _ in range(len(columns))]
        invalid_idxs += [idx]
    R += [r]

invalid_idxs = set(invalid_idxs)

cols_to_drop = joblib.load(os.path.join(checkpoints_dir, "cols_to_drop.pkl"))
imputer = joblib.load(os.path.join(checkpoints_dir,"imputer.pkl"))

R = pd.DataFrame(R, columns=columns)
R = R.drop(columns=cols_to_drop)
R = convert_to_float(R)
R = imputer.transform(R)

cols = [c for c in columns if c not in cols_to_drop]
cols = [c.lower() for c in cols]
cols = [c.replace("-", "_") for c in cols]
cols = [c.replace(" ", "_") for c in cols]
cols = [c.replace("(", "_") for c in cols]
cols = [c.replace(")", "") for c in cols]

cols_ = []
for c in cols:
    if c in cols_:
        cols_ += [c+"_1"]
    else:
        cols_ += [c]

print(len(cols_), "columns after processing")
print(len(set(cols_)), "unique columns after processing")

with open(outfile, "w") as f:
    writer = csv.writer(f)
    writer.writerow(cols_)
    for i, r in enumerate(R):
        if i in invalid_idxs:
            writer.writerow([None for _ in range(len(cols_))])
        else:
            writer.writerow(r)
