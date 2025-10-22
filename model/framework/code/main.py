import os
import csv
import json
import sys
import struct  # FIX: needed for read_smiles_bin
from rdkit import Chem
import pandas as pd
import numpy as np
from mordred import Calculator, descriptors
from timeout_decorator import timeout
from ersilia_pack_utils.core import write_out, read_smiles
import joblib

# ---------- CLI ----------
infile = sys.argv[1]
outfile = sys.argv[2]

# ---------- Paths / constants ----------
ROOT = os.path.abspath(os.path.dirname(__file__))
checkpoints_dir = os.path.join(ROOT, "..", "..", "checkpoints")
TIMEOUT_SEC = 60


# ---------- Descriptor calculator ----------
calc = Calculator(descriptors, ignore_3D=True)
# Reference list of all Mordred descriptor names
descriptor_columns = list(calc.pandas([Chem.MolFromSmiles("CCCC")]).columns)

@timeout(TIMEOUT_SEC)
def one_molecule(mol):
    # Return a flat list of descriptor values for this molecule
    # calc(mol) → Mordred result; convert to list to be safe
    return list(calc(mol))

# ---------- Main workflow ----------
# Read SMILES (ignore the input file's own columns)
_, smiles = read_smiles(infile)

rows = []
invalid_idxs = []

for idx, smi in enumerate(smiles):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError("RDKit failed to parse SMILES")
        row = one_molecule(mol)  # may timeout / raise
    except Exception:
        row = [np.nan] * len(descriptor_columns)  # use NaN, not None
        invalid_idxs.append(idx)
    rows.append(row)

# Build DataFrame with descriptor names
R = pd.DataFrame(rows, columns=descriptor_columns)

# Load preprocessing artifacts
cols_to_drop = joblib.load(os.path.join(checkpoints_dir, "cols_to_drop.pkl"))
imputer = joblib.load(os.path.join(checkpoints_dir, "imputer.pkl"))

# Drop columns safely (ignore missing in case lists differ)
R = R.drop(columns=[c for c in cols_to_drop if c in R.columns], errors="ignore")

# Coerce everything to numeric
R = R.apply(pd.to_numeric, errors="coerce")

# Impute → numpy array
R_imputed = imputer.transform(R)  # shape: (n_samples, n_features_after_drop)

# Clean output column names in the SAME order as R’s columns
clean_cols = []
seen = set()
for c in R.columns:
    c2 = (
        c.lower()
         .replace("-", "_")
         .replace(" ", "_")
         .replace("(", "_")
         .replace(")", "")
    )
    # ensure uniqueness
    if c2 in seen:
        suffix = 1
        while f"{c2}_{suffix}" in seen:
            suffix += 1
        c2 = f"{c2}_{suffix}"
    seen.add(c2)
    clean_cols.append(c2)

print(f"{len(clean_cols)} columns after processing")
print(f"{len(set(clean_cols))} unique columns after processing")
if invalid_idxs:
    print(f"Warning: {len(invalid_idxs)} invalid SMILES at rows {invalid_idxs[:10]}{'...' if len(invalid_idxs) > 10 else ''}")

# Write out (CSV or BIN), header = clean column names
write_out(R_imputed, clean_cols, outfile,'float32')
