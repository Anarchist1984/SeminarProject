import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs

def load_data(filepath):
    """Loads data from a CSV file."""
    return pd.read_csv(filepath)

def calculate_fingerprints(smiles_list):
    """Calculates fingerprints for a list of SMILES strings."""
    fingerprints = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        fp = Chem.RDKFingerprint(mol)
        fingerprints.append(fp)
    return fingerprints

def calculate_tanimoto_similarity(fp1, fp2):
    """Calculates Tanimoto similarity between two fingerprints."""
    return DataStructs.TanimotoSimilarity(fp1, fp2)
