import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pyopenms import MSExperiment, MzMLFile
from tqdm import tqdm


def select_top_peaks(mzs, intensities, num_peaks=5000):
    indices = np.argsort(intensities)[::-1][:num_peaks]
    return mzs[indices], intensities[indices]


def calculate_pairwise_distances(mzs, intensities):
    mzs = np.array(mzs)
    intensities = np.array(intensities)
    mz_diffs = np.abs(mzs[:, np.newaxis] - mzs)
    mz_diffs = np.triu(mz_diffs, k=1)
    mz_i, mz_j = np.where((mz_diffs > 0) & (mz_diffs <= 300))
    return list(
        zip(
            mzs[mz_i],
            mzs[mz_j],
            mz_diffs[mz_i, mz_j],
            intensities[mz_i] + intensities[mz_j],
        )
    )


def bin_and_sum_intensities(pairs):
    bins = np.arange(0, 300.005, 0.005)
    intensity_sums = np.zeros(len(bins))
    for mz1, mz2, diff, intensity_sum in pairs:
        idx = int(diff / 0.005)
        intensity_sums[idx] += intensity_sum
    return intensity_sums


def process_spectrum_batch(spectrum_batch):
    batch_results = []
    for mz_intensity_data in spectrum_batch:
        mzs, intensities = mz_intensity_data
        if len(mzs) > 1:
            mzs, intensities = select_top_peaks(mzs, intensities)
            pairs = calculate_pairwise_distances(mzs, intensities)
            batch_results.append(bin_and_sum_intensities(pairs))
        else:
            batch_results.append(np.zeros(60001))
    return np.sum(batch_results, axis=0)


def create_batches(data_list, batch_size):
    return [data_list[i : i + batch_size] for i in range(0, len(data_list), batch_size)]


def process_file(mzml_path):
    exp = MSExperiment()
    MzMLFile().load(mzml_path, exp)
    spectra = [
        (s.get_peaks()[0], s.get_peaks()[1])
        for s in exp.getSpectra()
        if s.getMSLevel() == 2
    ]
    batches = create_batches(spectra, 50)  # Batch size of 50 spectra
    with ProcessPoolExecutor() as executor:
        results = list(
            tqdm(
                executor.map(process_spectrum_batch, batches),
                total=len(batches),
                desc=f"Processing batches in {os.path.basename(mzml_path)}",
            )
        )
    total_intensities = np.sum(results, axis=0)
    return total_intensities


def save_to_csv(results, filename):
    bin_start = np.arange(0, 300.000, 0.005)
    bin_end = np.arange(0.005, 300.005, 0.005)

    # Ensure the results array does not exceed the bin array lengths
    results = results[: len(bin_start)]  # Truncate results if necessary

    # Create a DataFrame
    df = pd.DataFrame(
        {"Bin_Start": bin_start, "Bin_End": bin_end, "Summed_Intensities": results}
    )

    # Save to CSV
    df.to_csv(filename, index=False)
    print(f"Results saved to {filename}")


def load_from_csv(filename):
    df = pd.read_csv(filename)
    return df["Summed_Intensities"].values


if __name__ == "__main__":
    # Example usage
    mzml_dir = "PXD028735/"

    results_dict = {}
    files = [
        os.path.join(mzml_dir, f) for f in os.listdir(mzml_dir) if f.endswith(".mzML")
    ]

    for file in tqdm(files, desc="Processing Files"):
        if "ScanningSWATH" in file:
            continue

        csv_filename = os.path.join(
            mzml_dir, os.path.splitext(os.path.basename(file))[0] + "_intensities.csv"
        )
        if os.path.exists(csv_filename):
            print(f"Loading results from {csv_filename}")
            results_dict[os.path.basename(file)] = load_from_csv(csv_filename)
        else:
            result = process_file(file)
            results_dict[os.path.basename(file)] = result
            save_to_csv(result, csv_filename)

    # Print or process results stored in the dictionary
    for file, intensities in results_dict.items():
        print(f"Filename: {file}")
        print("Summed Intensities:", intensities)
        plt.scatter(np.arange(0, 300.000, 0.005), intensities, s=0.1, alpha=0.5)
        plt.savefig("img/" + file + ".png")
        plt.close()
