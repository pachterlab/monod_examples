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

Defining a Model
----------------------

To define a CME model of transcription and sequencing, initialize an instance of :py:class:`cme_toolbox.CMEModel`:

.. code-block:: python

 fitmodel = cme_toolbox.CMEModel('ProteinBursty','Poisson')

where ``biological_model = {'Bursty','Constitutive','Extrinsic','CIR'}`` represents the transcriptional process and ``sequencing_model = {'None','Bernoulli','Poisson'}`` represents the dynamics of the sampling process.


Inference
----------------

The entire inference procedure can be carried out in a single function,:py:func:`inference.perform_inference`, using anndata or an h5ad filename as input.

.. code-block:: python

 fitted_adata = inference.perform_inference(input_adata, fitmodel)

Additional keyword arguments can be used to adjust the settings of the fit. TODO: go through the rest of the keyword parameters.

To set the name of the output folder, set ``dataset_string='your_output_dirname'``. Set ``viz=True`` (default) to visualize the gene filtering.

You can optionally specifiy a gene length annotation filepath ``transcriptome_filepath``. To specify the number of genes to analyze, set ``n_genes``. 

A file ``gene_set.csv`` will be created, with the list of genes that meet filtering thresholds for all datasets, and the file ``genes.csv`` with the list of genes selected for further analysis. This gene list can be set manually using genes_to_fit.

If you do not want to use the default gene-filtering parameters for the specified model, you can manually set ``filt_param = {'min_means':min_means, 'max_maxes':max_maxes, 'min_maxes':min_maxes}``, where ``min_means`` is a list of the minimum allowed mean for selected genes, where each float in the list is for a different modality, with the modalities ordered as in ``CMEModel.model_modalities``. ``max_maxes`` and ``min_maxes`` are specified similarly, for the maximum and minimum maximum counts respectively for selected genes.

Monod will then iterate over all sampling parameter grid points using ``n_cores`` processors.

Fit Parameters
------------------

If you do not want to use the default fit parameters, you can specify them in the call to :py:func:`inference.perform_inference`, using the following keywords.

``phys_lb`` and ``phys_ub`` are bounds on the transcriptional process model parameters.
``samp_lb`` and ``samp_ub`` are bounds on the sampling process model parameters.
``gridsize`` defines the grid for the sampling parameter scan.
``gradient_params`` defines the gradient optimization parameters, e.g. ``gradient_params = {'max_iterations':5,'init_pattern':'moments','num_restarts':1}``

If a transcriptome length annotation is provided, lengths will be used in determining the nascent RNA capture rate (to model priming at ubiquitous internal polyA sites). If lengths are not given, the keyword argument: ``poisson_average_log_length`` specifies, in base 10, what the universal multiplier on the nascent capture rate should be.

Alternatively one can define the search parameters and cluster the data. This will run the `meK-Means <https://github.com/pachterlab/CGP_2023/>`_ clustering algorithm (see the paper `here <https://www.biorxiv.org/content/10.1101/2023.09.17.558131v2>`_). In this case, the user should give a value for ``mek_means_params'' = (``k``, ``epochs``), where ``k`` is the user-defined number of clusters to learn and ``epochs`` is the numbers of rounds to learn the data clusters. All other parameters remain the same. 


Post-processing and QC
----------------

.. code-block:: python

analysis.run_qc(fitted_adata)

To retrieve the SearchResults and SearchData object from the fitted anndata object, you can run:

.. code-block:: python

search_result, search_data = fitted_adata.uns['search_result'], fitted_adata.uns['search_data']

To identify the technical noise parameter optimum, call a method of a SearchResults object:

.. code-block:: python

 search_result.find_sampling_optimum()

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

TODO: Do we want to add noise decomposition?

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
 
