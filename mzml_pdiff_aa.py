import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pyopenms import MSExperiment, MzMLFile
from tqdm import tqdm

# List of delta m/z values for common amino acids (example values)
amino_acid_masses = {
    "G": 57.02146,
    "A": 71.03711,
    "S": 87.03203,
    "P": 97.05276,
    "V": 99.06841,
    "T": 101.04768,
    "C": 103.00919,
    "L": 113.08406,
    "I": 113.08406,
    "N": 114.04293,
    "D": 115.02694,
    "Q": 128.05858,
    "K": 128.09496,
    "E": 129.04259,
    "M": 131.04049,
    "H": 137.05891,
    "F": 147.06841,
    "R": 156.10111,
    "Y": 163.06333,
    "W": 186.07931,
    "C[Carbamidomethyl]": 57.0215,  # Example for "C with carbamidomethyl modification":
    "K[methyl]": 142.1267,  # Example for "K with methyl modification":
    "Q[Pyro-glu]": 128.0586,  # Example for "Q with pyro-glu modification":
    "K[Dimethyl]": 156.1406,  # Example for "K with dimethyl modification":
    "S[Phospho]": 166.9987,  # Example for "S with phospho modification":
    "T[Phospho]": 181.0147,  # Example for "T with phospho modification":
    "Y[Phospho]": 243.0297,  # Example for "Y with phospho modification":
    "K[Trimethyl]": 170.1545,  # Example for "K with trimethyl modification":
    "K[Acetyl]": 170.1055,  # Example for "K with acetyl modification":
    "N[Deamidated]": 115.0270,  # Example for "N with deamidated modification":
    "K[Propionyl]": 184.1194,  # Example for "K with propionyl modification":
    "K[Butyryl]": 198.1333,  # Example for "K with butyryl modification":
    "K[Formyl]": 168.0892,  # Example for "K with formyl modification":
    "K[Oxidation]": 144.1059,  # Example for "K with oxidation modification":
    "K[Malonyl]": 184.1055,  # Example for "K with malonyl modification":
    "K[Succinyl]": 198.1194,  # Example for "K with succinyl modification":
    "R[Methyl]": 170.0914,  # Example for "R with methyl modification":
    "R[Dimethyl]": 184.1053,  # Example for "R with dimethyl modification":
    "R[Trimethyl]": 198.1192,  # Example for "R with trimethyl modification":
    "Q[Deamidated]": 115.0270,  # Example for "Q with deamidated modification":
    "N[Oxidation]": 115.0269,  # Example for "N with oxidation modification":
    "M[Oxidation]": 147.0354,  # Example for "M with oxidation modification":
    "H[Oxidation]": 153.0192,  # Example for "H with oxidation modification":
    "M[Carbamylation]": 125.0477,  # Example for "M with carbamylation modification":
    "N[Amonia-loss]": 113.0473,  # Example for "N with ammonia-loss modification":
    "R[Amonia-loss]": 129.0426,  # Example for "R with ammonia-loss modification":
    "Q[Amonia-loss]": 114.0429,  # Example for "Q with ammonia-loss modification":
    "N[Pyro-glu]": 128.0586,  # Example for "N with pyro-glu modification":
    "K[Pyro-glu]": 128.0949,  # Example for "K with pyro-glu modification":
    "R[Pyro-glu]": 156.1011,  # Example for "R with pyro-glu modification":
}

"""

"""

"""

"""

# Create a new dictionary for two amino acid combinations
combined_masses = {}
for aa1 in amino_acid_masses:
    for aa2 in amino_acid_masses:
        combination_key = aa1 + aa2
        combination_mass = amino_acid_masses[aa1] + amino_acid_masses[aa2]
        combined_masses[combination_key] = combination_mass

"""
# Create a new dictionary for three amino acid combinations
triple_combined_masses = {}
for aa1 in amino_acid_masses:
    for aa2 in amino_acid_masses:
        for aa3 in amino_acid_masses:
            combination_key = aa1 + aa2 + aa3
            combination_mass = (
                amino_acid_masses[aa1] + amino_acid_masses[aa2] + amino_acid_masses[aa3]
            )
            triple_combined_masses[combination_key] = combination_mass
"""

# Combine the original amino acid masses with the new combination masses
all_masses = {**amino_acid_masses}  # , **combined_masses , **triple_combined_masses

# Converting to a sorted list of masses for easier indexing
sorted_masses = sorted(all_masses.items(), key=lambda x: x[1])
amino_acid_deltas = [mass for name, mass in sorted_masses]
amino_acid_deltas_names = [name for name, mass in sorted_masses]


def calculate_pairwise_distances(mzs, amino_acid_deltas, tolerance=0.005, max_mz=300.0):
    mzs = np.array(mzs)
    if len(mzs) < 2:  # Ensure there are at least two m/z values to compare
        return np.zeros(len(amino_acid_deltas), dtype=int)
    mz_diffs = np.abs(mzs[:, np.newaxis] - mzs)
    mz_diffs = np.triu(mz_diffs, k=1)
    mz_i, mz_j = np.where((mz_diffs > 0) & (mz_diffs <= max_mz))
    if not mz_i.size:  # No valid differences found
        return np.zeros(len(amino_acid_deltas), dtype=int)
    mz_diffs = mz_diffs[mz_i, mz_j]
    comparison = np.abs(mz_diffs[:, np.newaxis] - amino_acid_deltas) <= tolerance
    return np.sum(comparison, axis=0)


def bin_and_sum_intensities(pairs):
    bins = np.arange(0, 300.005, 0.005)
    intensity_sums = np.zeros(len(bins))
    for mz1, mz2, diff, intensity_sum in pairs:
        idx = int(diff / 0.005)
        intensity_sums[idx] += intensity_sum
    return intensity_sums


def process_spectrum_batch(spectrum_batch, amino_acid_deltas):
    results_matrix = []
    for mzs, intensities in spectrum_batch:
        if len(mzs) > 1:
            results_matrix.append(calculate_pairwise_distances(mzs, amino_acid_deltas))
        else:
            results_matrix.append(
                np.zeros(len(amino_acid_deltas), dtype=int)
            )  # Ensure consistent shape
    return np.array(results_matrix)


def create_batches(data_list, batch_size):
    return [data_list[i : i + batch_size] for i in range(0, len(data_list), batch_size)]


def process_file(mzml_path, amino_acid_deltas):
    exp = MSExperiment()
    MzMLFile().load(mzml_path, exp)
    spectra = [
        (s.get_peaks()[0], s.get_peaks()[1])  # select_top_peaks
        for s in exp.getSpectra()
        if s.getMSLevel() == 2
    ]

    batches = create_batches(spectra, 100)  # Assuming batch size of 500

    with ProcessPoolExecutor() as executor:
        matrix_list = list(
            tqdm(
                executor.map(
                    process_batch_wrapper, batches, [amino_acid_deltas] * len(batches)
                ),
                total=len(batches),
                desc="Processing batches",
            )
        )
        matrix_list = [
            m for m in matrix_list if m.size != 0
        ]  # Filter out empty matrices

    if matrix_list:  # Check if there's anything to concatenate
        final_matrix = np.vstack(matrix_list)
    else:
        final_matrix = np.array([])  # Or handle empty results appropriately

    return final_matrix


def save_to_csv(data, filename):
    # Convert the numpy array to a pandas DataFrame
    # df = pd.DataFrame(data)

    # Optional: If you want to add column headers
    # df.columns = ['Column1', 'Column2', 'Column3', ..., 'Column20']  # Adjust column names as necessary

    # Write the DataFrame to a CSV file
    data.to_csv(
        filename, index=False
    )  # Set index=False if you do not want to write row numbers


def load_from_csv(filename):
    df = pd.read_csv(filename)
    return df["Summed_Intensities"].values


def process_batch_wrapper(batch, amino_acid_deltas):
    return process_spectrum_batch(batch, amino_acid_deltas)


if __name__ == "__main__":
    # Example usage
    mzml_dir = "PXD028735/"

    results_dict = {}
    files = [
        os.path.join(mzml_dir, f)
        for f in os.listdir(mzml_dir)
        if f.endswith(".mzML")  # mzML
    ]

    to_analyze = [
        "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01.mzML",
        "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.mzML",
        "LFQ_Orbitrap_AIF_Human_01.mzML",
        "LFQ_Orbitrap_DDA_Human_01.mzML",
    ]

    for file in tqdm(files, desc="Processing Files"):
        # if "part_3939.327224524804_4118.387552912295.mzml" not in file:
        #    continue
        if "SWATH" in file:
            continue
        if "ScanningSWATH" in file:
            continue
        if "_GP_" in file:
            continue
        if file.split("/")[-1] not in to_analyze:
            continue
        # if "Orbitrap" not in file:
        #    continue

        print(file)

        csv_filename = os.path.join(
            mzml_dir,
            os.path.splitext(os.path.basename(file))[0]
            + "_intensities_aa_delta_single_mod.csv",  # double triple
        )
        if os.path.exists(csv_filename) and False:
            print(f"Loading results from {csv_filename}")
            # results_dict[os.path.basename(file)] = load_from_csv(csv_filename)
        else:
            result = process_file(file, amino_acid_deltas)
            result = pd.DataFrame(result)
            result.columns = amino_acid_deltas_names
            results_dict[os.path.basename(file)] = result
            save_to_csv(result, csv_filename)

    # Print or process results stored in the dictionary
    for file, intensities in results_dict.items():
        print(f"Filename: {file}")
        print("Summed Intensities:", intensities)
        plt.scatter(np.arange(0, 300.000, 0.005), intensities, s=0.1, alpha=0.5)
        plt.savefig(file + ".png")
        plt.close()
