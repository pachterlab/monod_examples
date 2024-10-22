Usage
=====

.. _installation:

Installation
------------

To use *Monod*, install it from `pip` (TODO make pip-installable):

.. code-block:: console

 pip install monod
 
To use it in your code, import the package components:

.. code-block:: python

 import monod
 from monod import *
 from monod.extract_data import *
 from monod.cme_toolbox import CMEModel
 from monod import inference, mminference
 from monod.analysis import *

Quantification 
----------------
To generate an intronic index for ``kb``, obtain a reference genome and run ``kb ref``:

.. code-block:: console

 kb ref -i ./refdata-gex-mm10-2020-A/kallisto/index.idx 
 -g ./refdata-gex-mm10-2020-A/t2g_mm10.txt 
 -f1 ./refdata-gex-mm10-2020-A/kallisto/cdna.fa 
 -f2 ./refdata-gex-mm10-2020-A/kallisto/intron.fa 
 -c1 ./refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt 
 -c2 ./refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt 
 --workflow lamanno 
 ./refdata-gex-mm10-2020-A/fasta/genome.fa 
 ./refdata-gex-mm10-2020-A/genes/genes.gtf
 
To then generate the count matrices in ``h5ad`` format, use ``kb count``:

.. code-block:: console

   kb count --verbose 
   -i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx 
   -g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt 
   -x 10xv3 
   -o OUTDIR/ 
   -t 30 -m 30G 
   -c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt 
   -c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt 
   --workflow lamanno --filter bustools --overwrite --h5ad 
   DATASET_FASTQ_LOCATIONS

This generates an ``h5ad`` file in ``OUTDIR/counts_filtered/``.

Pre-processing 
----------------

Additional cell/gene filtering can be done using scanpy, and the resulting anndata object can be directly input into Monod. Alternatively, the raw filename can be given.

TODO: add scanpy code.

Inference
----------------

The entire inference procedure can be carried out in a single function,:py:func:`inference.perform_inference`, using anndata or an h5ad filename as input.

.. code-block:: python

fitted_adata = inference.perform_inference(input_adata, fitmodel)

Additional keyword arguments can be used to adjust the settings of the fit. TODO: go through the rest of the keyword parameters.

.. code-block:: python

, n_genes=n_genes, seed=5,
            phys_lb=lb, phys_ub=ub, gridsize=grid, samp_lb=samp_lb, samp_ub=samp_ub,
            gradient_params = {'max_iterations':5,'init_pattern':'moments','num_restarts':1},
                                         transcriptome_filepath=transcriptome_filepath, poisson_average_log_length=5, dataset_string='cite_fit', viz=True,
                                         num_cores=32, mek_means_params=mek_means_params, filt_param = {'min_means':[0, 0, 0], 'max_maxes':[3500, 3500, 1000000], 'min_maxes':[0,0, 0]}

To set the name of the output folder, set dataset_name='your_output_dirname'.

You can optionally specifiy a gene length annotation filepath ``transcriptome_filepath``. To specify the number of genes to analyze, set ``n_genes``. 

A file ``gene_set.csv`` will be created, with the list of genes that meet filtering thresholds for all datasets, and the file ``genes.csv`` with the list of genes selected for further analysis. This gene list can be set manually using genes_to_fit.

Model, data, and parameter definition 
----------------

To define a CME model of transcription and sequencing, initialize an instance of :py:class:`cme_toolbox.CMEModel`:

.. code-block:: python

 CMEModel(biological_model,sequencing_model)

where ``biological_model = {'Bursty','Constitutive','Extrinsic','CIR'}`` represents the transcriptional process and ``sequencing_model = {'None','Bernoulli','Poisson'}`` represents the dynamics of the sampling process.

To define the search parameters, initialize an instance of :py:class:`inference.InferenceParameters`:

.. code-block:: python

 inference_parameters = inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                     dataset_string,fitmodel,use_lengths)

where ``phys_lb`` and ``phys_ub`` are bounds on the transcriptional process model parameters, ``samp_lb`` and ``samp_ub`` are bounds on the sampling process model parameters, ``gridsize`` defines the grid for the sampling parameter scan, and ``use_lengths`` determines whether the unspliced mRNA capture rate depends on the gene length (to model priming at ubiquitous internal polyA sites).

Alternatively one can define the search parameters and cluster the data. This will run the `meK-Means <https://github.com/pachterlab/CGP_2023/>`_ clustering algorithm (see the paper `here <https://www.biorxiv.org/content/10.1101/2023.09.17.558131v2>`_). We initialize an instance of :py:class:`mminference.InferenceParameters`:

.. code-block:: python

 inference_parameters = mminference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                     k,epochs,\
                     dataset_string,fitmodel,use_lengths)

where ``k`` is the user-defined number of clusters to learn and ``epochs`` is the numbers of rounds to learn the data clusters. All other parameters remain the same as :py:class:`inference.InferenceParameters`. 

To create a ``SearchData`` object to input into the inference process, run :py:func:`extract_data.extract_data`:

.. code-block:: python

 search_data = extract_data(loom_filepath, transcriptome_filepath, dataset_name,
                                dataset_string, dir_string)

Running the inference pipeline 
----------------

To run the pipeline, simply call the following parallelized code:

.. code-block:: python

  result_string = inference_parameters.fit_all_grid_points(n_cores,search_data)

This will iterate over all grid points using ``n_cores`` processors.

Post-processing and QC
----------------
To load the search results, import the file string:

.. code-block:: python

 sr = load_search_results(result_string)

To identify the technical noise parameter optimum, call a method of a SearchResults object `sr`:

.. code-block:: python

 sr.find_sampling_optimum()

Optionally, test its stability under subsampling and chi-squared testing:

.. code-block:: python

 fig1,ax1 = plt.subplots(1,1)
 sr.plot_landscape(ax1)
 _=sr.chisquare_testing(sd)
 sr.resample_opt_viz()
 sr.resample_opt_mc_viz()
 sr.chisq_best_param_correction(sd,viz=True)

Optionally, examine whether the distribution fits match the raw data:

.. code-block:: python

 sr.plot_gene_distributions(sd,marg='joint')
 sr.plot_gene_distributions(sd,marg='nascent')
 sr.plot_gene_distributions(sd,marg='mature')

To chracterize the uncertainty, variation, and bias in biological parameters, compute the standard errors of their maximum likelihood estimates, then plot their distributions and dependence on length (which should be minimal):

.. code-block:: python

 sr.compute_sigma(sd,num_cores)
 sr.plot_param_L_dep(plot_errorbars=True,plot_fit=True)
 sr.plot_param_marg()

As the standard error computation is typically computationally intensive, it is useful to store an updated copy on disk after evaluating it:

.. code-block:: python

 sr.update_on_disk()

Noise decomposition
----------------
Two complementary methods are available for investigating the contributions of different noise sources. The first is non-parametric; calling a method of a ``SearchData`` object ``sd`` returns the fractions of variation retained and discarded after normalization and log-transformation:

.. code-block:: python

 f = sd.get_noise_decomp()
 
These fractions are not guaranteed to be positive, because the transformations may *increase* the relative spread of the data. On the other hand, if a fit has been completed, a method of a ``SearchResults`` object ``sr`` reports the fractions of intrinsic, extrinsic (bursty), and technical variation for each gene and molecular species:

.. code-block:: python

 f = sr.get_noise_decomp()
 
Differential parameter value identification
----------------
Given a set of matched datasets, run with the same model over the same set of genes, two approaches are available for identifying putative patterns of differential expression and regulation. A moment-based, biology-agnostic one uses a simple *t*-test to identify differences in the means of spliced counts in ``SearchData`` objects ``sd1`` and ``sd2``:

.. code-block:: python

 gf = compute_diffexp(sd1,sd2)
 
where ``gf`` is boolean vector that reports ``True`` if the gene is identified as DE. However, this approach cannot identify differences if biological parameters change in a correlated way and the mean stays the same. We introduce a more mechanistic approach for the identification of differential expression suggested by parameter variation, based on two ``SearchResults`` objects ``sr1`` and ``sr2``:

.. code-block:: python

 gf = compute_diffreg(sr1,sr2)
 
where ``gf`` is a two-dimensional boolean array that reports ``True`` if a particular *parameter* is identified as DE. After using these arrays to find a subpopulation of interest -- e.g., genes that do not exhibit variation in the spliced mean, but do exhibit modulation in the burst size -- it is possible to plug the gene filter ``genes_to_plot`` back in to inspect the raw data and fits:

.. code-block:: python

 gf = compare_gene_distributions(sr_arr,sd_arr,genes_to_plot=genes_to_plot)
 
