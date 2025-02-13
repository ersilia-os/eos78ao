import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
import joblib

def remove_cols(df):
    columns_to_drop = []
    for col in df.columns:
        num_total = len(df[col])  # Total number of rows in column
        num_numeric = pd.to_numeric(df[col], errors='coerce').notna().sum()  # Count valid numeric values
        # If more than `threshold` proportion of values are strings, mark for removal
        if num_numeric / num_total < (1 - 0.5):  
            print(f"Column '{col}' has mostly string values and will be dropped.")
            columns_to_drop.append(col)
    return columns_to_drop

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

def clean_up(df):
    # Drop molecules where more than 20% of values are NaN
    threshold = 0.2 * df.shape[1]
    df_cleaned = df.dropna(thresh=df.shape[1] - threshold)
    return df_cleaned

def impute(df):    
    imputer = SimpleImputer(strategy='median')
    imputer.fit(df)
    return imputer

df = pd.read_csv("drugbank_mordred.csv")
print(df.shape)
cols_to_drop = remove_cols(df)
joblib.dump(cols_to_drop, "cols_to_drop.pkl")
print(len(cols_to_drop))
df = df.drop(columns=cols_to_drop)
print(df.shape)
df = convert_to_float(df)
df = clean_up(df)
print(df.shape)

imputer = impute(df)
joblib.dump(imputer, 'imputer.pkl')
df_imputed = imputer.transform(df)
df_imputed = pd.DataFrame(df_imputed, columns=df.columns)
df_imputed.to_csv("drugbank_mordred_imputed.csv", index=False)
