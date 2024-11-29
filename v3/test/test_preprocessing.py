import unittest
import pandas as pd
from src.preprocessing import load_data, calculate_fingerprints, calculate_tanimoto_similarity
from rdkit import Chem

class TestPreprocessing(unittest.TestCase):
    def test_load_data(self):
        # Create a sample CSV file
        data = {'col1': [1, 2], 'col2': [3, 4]}
        df = pd.DataFrame(data)
        df.to_csv('test.csv', index=False)

        # Load the data
        loaded_data = load_data('test.csv')

        # Assert that the data is loaded correctly
        self.assertEqual(len(loaded_data), 2)
        self.assertEqual(len(loaded_data.columns), 2)

    def test_calculate_fingerprints(self):
        smiles_list = ['CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C']
        fingerprints = calculate_fingerprints(smiles_list)
        self.assertEqual(len(fingerprints), 2)
        for fp in fingerprints:
            self.assertTrue(isinstance(fp, list))

    def test_calculate_tanimoto_similarity(self):
        mol1 = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
        mol2 = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
        similarity = calculate_tanimoto_similarity(fp1, fp2)
        self.assertTrue(0 <= similarity <= 1)

if __name__ == '__main__':
    unittest.main()
