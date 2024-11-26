#@title Library Function
import os
import math
import pandas as pd
from rdkit.Chem import PandasTools
from chembl_webresource_client.new_client import new_client


class UniProtDataProcessor:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        self.targets_api = new_client.target
        self.compounds_api = new_client.molecule
        self.bioactivities_api = new_client.activity

    def fetch_target_chembl_id(self):
        """Fetches the ChEMBL ID for the given UniProt ID."""
        targets = self.targets_api.get(
            target_components__accession=self.uniprot_id
        ).only("target_chembl_id")
        targets_df = pd.DataFrame.from_records(targets)
        if targets_df.empty:
            raise ValueError("No targets found for the given UniProt ID.")
        return targets_df.iloc[0]["target_chembl_id"]

    def fetch_bioactivities(self, chembl_id):
        """Fetches bioactivities for the target ChEMBL ID."""
        bioactivities = self.bioactivities_api.filter(
            target_chembl_id=chembl_id, type="IC50", relation="=", assay_type="B"
        ).only(
            "molecule_chembl_id",
            "standard_value",
            "standard_units",
        )
        bioactivities_df = pd.DataFrame.from_dict(bioactivities)
        if bioactivities_df.empty:
            raise ValueError("No bioactivities found for the target ChEMBL ID.")
        bioactivities_df = bioactivities_df.astype({"standard_value": "float64"})
        bioactivities_df.dropna(inplace=True)
        bioactivities_df = bioactivities_df[bioactivities_df["standard_units"] == "nM"]
        bioactivities_df.drop_duplicates("molecule_chembl_id", inplace=True)
        bioactivities_df.rename(columns={"standard_value": "IC50"}, inplace=True)
        return bioactivities_df

    def fetch_compounds(self, molecule_ids):
        """Fetches compound information for given molecule ChEMBL IDs."""
        compounds_provider = compounds_api.filter(
        molecule_chembl_id__in=list(bioactivities_df["molecule_chembl_id"])
            ).only("molecule_chembl_id", "molecule_structures")
        compounds = list(tqdm(compounds_provider))
        compounds_df = pd.DataFrame.from_records(compounds)
        compounds_df.dropna(inplace=True)
        compounds_df["smiles"] = compounds_df["molecule_structures"].apply(
            lambda x: x.get("canonical_smiles") if x else None
        )
        compounds_df.dropna(subset=["smiles"], inplace=True)
        compounds_df.drop("molecule_structures", axis=1, inplace=True)
        return compounds_df

    def calculate_pic50(self, df):
        """Calculates pIC50 values and adds them to the DataFrame."""
        def convert_ic50_to_pic50(ic50_value):
            return 9 - math.log10(ic50_value) if ic50_value > 0 else float("nan")

        df["pIC50"] = df["IC50"].apply(convert_ic50_to_pic50)
        return df

    def process_data(self):
        """Main method to fetch and process data for the UniProt ID."""
        chembl_id = self.fetch_target_chembl_id()
        print(f"Target ChEMBL ID: {chembl_id}")
        bioactivities_df = self.fetch_bioactivities(chembl_id)
        print(f"Bioactivities shape: {bioactivities_df.shape}")
        compounds_df = self.fetch_compounds(bioactivities_df["molecule_chembl_id"])
        print(f"Compounds shape: {compounds_df.shape}")
        output_df = pd.merge(
            bioactivities_df, compounds_df, on="molecule_chembl_id"
        )
        output_df = self.calculate_pic50(output_df)
        PandasTools.AddMoleculeColumnToFrame(output_df, smilesCol="smiles")
        output_df = output_df.drop("ROMol", axis=1)
        print(f"Final output shape: {output_df.shape}")
        output_df.sort_values(by="pIC50", ascending=False, inplace=True)
        output_df.reset_index(drop=True, inplace=True)
        return output_df


# Example usage
uniprot_id = "P14780"
processor = UniProtDataProcessor(uniprot_id)
processed_data = processor.process_data()
processed_data.head()