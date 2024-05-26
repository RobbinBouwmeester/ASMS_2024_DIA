# ASMS_2024_DIA

All the files for analysis are prepared from the raw files. Only if you want to follow the analysis from the start to the end first download the following files (repository: [pride/PXD028735](https://www.ebi.ac.uk/pride/archive/projects/PXD028735)):

* LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01.raw
* LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.raw
* LFQ_Orbitrap_AIF_Human_01.raw
* LFQ_Orbitrap_DDA_Human_01.raw
* LFQ_timsTOFPro_diaPASEF_Human_01.d
* LFQ_timsTOFPro_PASEF_Human_01.d
* LFQ_TTOF6600_SWATH_Human_01.wiff.scan
* LFQ_TTOF6600_SWATH_Human_01.wiff
* LFQ_TTOF6600_DDA_Human_01.wiff.scan
* LFQ_TTOF6600_DDA_Human_01.wiff

Conver the files to mzml using [MSConvert](https://proteowizard.sourceforge.io/download.html) with vendor peak picking and put them in the mzml file folder. For the timsTOF files please use [tdf2mzml](https://github.com/mafreitas/tdf2mzml).

# Analysis
The following section will explain per notebook or python script what analysis is done and what the resulting figures are.

## mzml_pdiff_aa.py
Calculate delta m/z distances for all the MS2 spectra. This is later used to quantify potential ambiguity in DDA VS DIA.

## diann_searches_mods.ipynb
Do a comparison of different DIA-NN searches and the peptidoforms it identified. The variable modifications are changed between the searches.

## FragmentCollision.ipynb
Calculate possible isobaric overlapping amino acid and modification combinations for fragment peaks.

## read_speclib_diann.ipynb
Read the spectral library from DIA-NN to spot potential differences for different peptidoforms.

## search_spacesize.ipynb
Calculate the search space size based on set parameters and fasta.

## visualize_aa_diff.ipynb
Visualize the delta mz values that correspond to amino acids and modifications.

## visualize_ms1.ipynb
Visualize the number of MS1 peaks in DDA and DIA data.

## visualize_ms2.ipynb
Visualize the number of MS2 peaks in DDA and DIA data.

## histone_analysis/spectral_library_calc.ipynb
Compare predictions from DIA-NN to observed histone peptide fragment intensity spectra.
