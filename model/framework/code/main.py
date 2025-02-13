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

TIMEOUT_SEC=3

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
                    #print(f"Non-numeric value found at row {index}, column '{col}, {value}'")
                    df.at[index, col] = np.nan
    df = df.applymap(lambda x: np.nan if pd.isna(x) else x)
    return df


#columns = calc._name_dict.keys()
columns = list(calc.pandas([Chem.MolFromSmiles("CCCC")]).columns)

R = []
for smi in smiles:
    try:
        r = one_molecule(Chem.MolFromSmiles(smi))
    except:
        r = [None for _ in range(len(columns))]
    R += [r]

cols_to_drop = joblib.load(os.path.join(checkpoints_dir, "cols_to_drop.pkl"))
imputer = joblib.load(os.path.join(checkpoints_dir,"imputer.pkl"))

R = pd.DataFrame(R, columns=columns)
R = R.drop(columns=cols_to_drop)
R = convert_to_float(R)
R = imputer.transform(R)

cols = [c for c in columns if c not in cols_to_drop]

#TODO: we are now imputing molecules that might be all NaNs originally

with open(outfile, "w") as f:
    writer = csv.writer(f)
    writer.writerow(cols)
    for r in R:
        writer.writerow(r)
