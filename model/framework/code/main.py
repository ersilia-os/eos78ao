import os
import csv
import json
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
TIMEOUT_SEC=60

def read_smiles_csv(in_file):
  with open(in_file, "r") as f:
    reader = csv.reader(f)
    cols = next(reader)
    data = [r[0] for r in reader]
    return cols, data

def read_smiles_bin(in_file):
  with open(in_file, "rb") as f:
    data = f.read()

  mv = memoryview(data)
  nl = mv.tobytes().find(b"\n")
  meta = json.loads(mv[:nl].tobytes().decode("utf-8"))
  cols = meta.get("columns", [])
  count = meta.get("count", 0)

  smiles_list = [None] * count
  offset = nl + 1
  for i in range(count):
    (length,) = struct.unpack_from(">I", mv, offset)
    offset += 4
    smiles_list[i] = mv[offset : offset + length].tobytes().decode("utf-8")
    offset += length

  return cols, smiles_list
    
def read_smiles(in_file):
  if in_file.endswith(".bin"):
    return read_smiles_bin(in_file)
  return read_smiles_csv(in_file)

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

def write_out_csv(results, header, file):
  with open(file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for r in results:
      writer.writerow(r)


def write_out_bin(results, header, file):
  arr = np.asarray(results, dtype=np.float32)
  meta = {"columns": header, "shape": arr.shape, "dtype": "float32"}
  meta_bytes = (json.dumps(meta) + "\n").encode("utf-8")

  with open(file, "wb") as f:
    f.write(meta_bytes)
    f.truncate(len(meta_bytes) + arr.nbytes)

  m = np.memmap(
    file, dtype=arr.dtype, mode="r+", offset=len(meta_bytes), shape=arr.shape
  )
  m[:] = arr
  m.flush()

def write_out(results, header, file):
  if file.endswith(".bin"):
    write_out_bin(results, header, file)
  elif file.endswith(".csv"):
    write_out_csv(results, header, file)
  else:
    raise ValueError(f"Unsupported extension for {file!r}")

calc = Calculator(descriptors, ignore_3D=True)

columns = list(calc.pandas([Chem.MolFromSmiles("CCCC")]).columns)

smiles = read_smiles(infile)

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

write_out(R, cols_, outfile)
