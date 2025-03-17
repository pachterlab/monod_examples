This directory contains notebooks used to generate manuscript figures and results. 

* `Fig1`: Folder containing scripts to generate subplots for Figure 1.
  * `gg230313_desai_fit.ipynb` to fit Desai data.
  * `gg230328_diffexp.ipynb` to summarize Figure 4 of "Length Biases..." and the Desai fits.
  * `tc230303_germcell_fit_analysis.ipynb` to fit and analyze germ cell data. 
  * `fetchGermCells.py` and `concatGermCells.py` to fetch and organize germ cell data.

* `Fig2`: Folder containing scripts to generate subplots for Figure 2.
  * `tc20241021_HCA_cancer_fit_analysis.ipynb` to fit and analyze HCA cancer data.

* `Fig3`: Folder containing scripts to generate subplots for Figure 3.
  * `process_Lu.py` to process radiation treatment data before fitting. 
  * `mc241031_fitLu.py` to fit radiation treatment data.
  * `analysisLu.ipynb` to analyze radiation treatment data.
  *  `memoryProfiling_Lu_T.py`and `memoryProfiling_Lu_T.sh` run Monod several times with different number of cells to test memory requirements. 

* `Fig4`: Folder containing scripts to generate subplots for Figure 4.
  * `mc240513_allen_B08_fit.py` to  fit Allen sample B08 cell types with different models to compare biophysical models.
  * `modelComparisonAllen.py` to extract metrics of model comparison from fits. 
  * `modelComparisonsAllen_analysis.ipynb' to further compare models and generate figures. 


* `Fig5`: Folder containing scripts to generate subplots for Figure 5.
  * `compareKhatebTakeiBurstsizes.ipynb` to compare the Khateb fit and Takei reported burstsizes per gene. 
  * `processKhateb.ipynb` process data from mouse ESCs before fitting. 
  * `fitKhatebMonod.py` fit data from mouse ESCs. 
  * `processTakei.ipynb` process fluorescence measurement data from mouse ESCs before fitting. 
  * `gg230316_brain_nuc_fit.ipynb` to fit the single-cell and single-nucleus data from 10x.
  * `gg230316_brain_nuc_analysis.ipynb` to analyze the results. 


* `Supplemental`: Folder containing scripts to generate tables and figures for the supplement. 
  * `gg220623_moments.m` is a MATLAB script used to calculate lower moments of the 18 implemented models (6 of biology, 3 of technical noise).
  * `gg230327_allen_cellsubtypes_fit.ipynb` to fit the glutamatergic subtypes and cell type using Monod.
  * `gg230328_allen_cellsubtypes_norm_analysis.ipynb` to analyze the performance of normalization and dimensionality reduction techniques, and compare them to the Monod results.
 
	
The raw data and full search results have been deposited on Zenodo. The prerequisite data and fits corresponding to "Length Biases..."  content is available at [Zenodo](https://zenodo.org/record/7388133) and [the manuscript GitHub](https://github.com/pachterlab/GP_2021_3).

The `figs/` directory contains manuscript figures.

`celltype_splitter.ipynb` extracts glutamatergic and GABAergic cell types for `sample_data`.
