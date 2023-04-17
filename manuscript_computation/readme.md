This directory contains notebooks used to generate manuscript figures and results. 

* Figures 1 and 2:
  * `gg230327_allen_cellsubtypes_fit.ipynb` to fit the glutamatergic subtypes and cell type using *Monod*.
  * `gg230328_allen_cellsubtypes_norm_analysis.ipynb` to analyze the performance of normalization and dimensionality reduction techniques, and compare them to the *Monod* results.
* Figure 3:
  * `gg230316_brain_nuc_fit.ipynb` to fit the single-cell and single-nucleus data from 10x.
  * `gg230328_brain_nuc_analysis.ipynb` to analyze the results.
* Figure 4:
  * `gg230313_desai_fit.ipynb` to fit the Desai data.
  * `gg230328_diffexp.ipynb` to summarize Figure 4 of "Length Biases..." and the Desai fits.

The raw data and full search results have been deposited on Zenodo. The prerequisite data and fits corresponding to "Length Biases..."  content is available at [Zenodo](https://zenodo.org/record/7388133) and [the manuscript GitHub](https://github.com/pachterlab/GP_2021_3).

In addition, `gg220623_moments.m` is a MATLAB script used to calculate lower moments of the 18 implemented models (6 of biology, 3 of technical noise).

The `figs/` directory contains manuscript figures.

`celltype_splitter.ipynb` extracts glutamatergic and GABAergic cell types for the `sample_data`.
