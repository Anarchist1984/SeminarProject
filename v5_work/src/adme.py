import pandas as pd
import math
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors

class ADMETAnalyzer:
    def __init__(self, thresholds, scaled_threshold, properties_labels, y_max=None):
        self.thresholds = thresholds
        self.scaled_threshold = scaled_threshold
        self.properties_labels = properties_labels
        self.y_max = y_max

    def filter_ro5(self, molecules):
        """
        Filter molecules based on Lipinski's Rule of Five (Ro5).

        Parameters
        ----------
        molecules : pd.DataFrame
            DataFrame with a "smiles" column containing molecular SMILES strings.

        Returns
        -------
        tuple
            Two DataFrames: one with Ro5-compliant molecules and one with non-compliant molecules.
        """
        ro5_properties = molecules["smiles"].apply(self.calculate_ro5_properties)
        molecules = pd.concat([molecules, ro5_properties], axis=1)
        molecules_ro5_fulfilled = molecules[molecules["ro5_fulfilled"]]
        molecules_ro5_violated = molecules[~molecules["ro5_fulfilled"]]
        return molecules_ro5_fulfilled, molecules_ro5_violated

    def visualize_admet(self, dataframe):
        """
        Visualize ADMET properties of molecules using a radar chart.

        Parameters
        ----------
        dataframe : pd.DataFrame
            DataFrame containing molecular properties (molecular weight, n_hba, n_hbd, logp).
        """
        stats = self.calculate_mean_std(dataframe[["molecular_weight", "n_hba", "n_hbd", "logp"]])
        self.plot_radar(stats)

    @staticmethod
    def calculate_ro5_properties(smiles):
        """
        Calculate properties relevant to Lipinski's Rule of Five (Ro5) for a molecule.

        Parameters
        ----------
        smiles : str
            SMILES string of the molecule.

        Returns
        -------
        pandas.Series
            Molecular properties and Ro5 compliance.
        """
        molecule = Chem.MolFromSmiles(smiles)
        molecular_weight = Descriptors.ExactMolWt(molecule)
        n_hba = Descriptors.NumHAcceptors(molecule)
        n_hbd = Descriptors.NumHDonors(molecule)
        logp = Descriptors.MolLogP(molecule)
        conditions = [molecular_weight <= 500, n_hba <= 10, n_hbd <= 5, logp <= 5]
        ro5_fulfilled = sum(conditions) >= 3
        return pd.Series(
            [molecular_weight, n_hba, n_hbd, logp, ro5_fulfilled],
            index=["molecular_weight", "n_hba", "n_hbd", "logp", "ro5_fulfilled"],
        )

    @staticmethod
    def calculate_mean_std(dataframe):
        """
        Calculate mean and standard deviation for molecular properties.

        Parameters
        ----------
        dataframe : pd.DataFrame
            DataFrame containing molecular properties.

        Returns
        -------
        pd.DataFrame
            DataFrame with mean and standard deviation of properties.
        """
        return dataframe.describe().T[["mean", "std"]]

    @staticmethod
    def _scale_by_thresholds(stats, thresholds, scaled_threshold):
        for property_name in stats.index:
            if property_name not in thresholds:
                raise KeyError(f"Property '{property_name}' is missing from thresholds.")
        return stats.apply(lambda x: x / thresholds[x.name] * scaled_threshold, axis=1)

    @staticmethod
    def _define_radial_axes_angles(n_axes):
        return [i / float(n_axes) * 2 * math.pi for i in range(n_axes)] + [0]

    def plot_radar(self, stats):
        """
        Plot a radar chart of molecular properties.

        Parameters
        ----------
        stats : pd.DataFrame
            DataFrame with mean and standard deviation of molecular properties.
        """
        x_angles = self._define_radial_axes_angles(len(stats))
        stats_scaled = self._scale_by_thresholds(stats, self.thresholds, self.scaled_threshold)
        stats_scaled = pd.concat([stats_scaled, stats_scaled.head(1)])

        plt.figure(figsize=(6, 6))
        ax = plt.subplot(111, polar=True)

        ax.fill(x_angles, [self.scaled_threshold] * len(x_angles), "cornflowerblue", alpha=0.2)
        ax.plot(x_angles, stats_scaled["mean"], "b", lw=3, ls="-")
        ax.plot(x_angles, stats_scaled["mean"] + stats_scaled["std"], "orange", lw=2, ls="--")
        ax.plot(x_angles, stats_scaled["mean"] - stats_scaled["std"], "orange", lw=2, ls="-.")

        ax.set_theta_offset(math.pi / 2)
        ax.set_theta_direction(-1)
        plt.xticks(x_angles[:-1], [])
        y_max = self.y_max or int(ax.get_yticks()[-1])
        plt.ylim(0, y_max)

        for angle, label in zip(x_angles[:-1], self.properties_labels):
            ha = "center" if angle in (0, math.pi) else ("left" if angle < math.pi else "right")
            ax.text(angle, y_max + 1, label, size=16, ha=ha)

        ax.legend(
            ("mean", "mean + std", "mean - std", "rule of five area"),
            loc=(1.1, 0.7), fontsize=16,
        )
        plt.show()

# Example usage:
# Define thresholds and labels
thresholds = {"molecular_weight": 500, "n_hba": 10, "n_hbd": 5, "logp": 5}
scaled_threshold = 5
properties_labels = ["Molecular weight (Da)", "# HBA", "# HBD", "LogP"]
y_max = 8

# Initialize the analyzer
analyzer = ADMETAnalyzer(thresholds, scaled_threshold, properties_labels, y_max)

# Example dataset (replace 'folder_path' with your dataset path)
molecules = pd.read_csv("folder_path.csv")
fulfilled, violated = analyzer.filter_ro5(molecules)
analyzer.visualize_admet(fulfilled)
