Advanced Options
=====

.. _fitparameters:

Additional Fit Parameters
------------

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

Analysis
----------------

To perform specific analyses on the fits, we can return to the Monod specific SearchResults and SearchData objects. To retrieve these from the fitted anndata object, you can run:

.. code-block:: python

  search_result, search_data = fitted_adata.uns['search_result'], fitted_adata.uns['search_data']

For fits with meK-Means clustering, we also have to specify a cluster for the search_result:

  search_result = monod_adata.uns['search_result_list'][cluster]

TODO: check the following are still correct.

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

Differential Parameters
----------------

Some more thoughts on differential parameters:

Given a set of matched datasets, run with the same model over the same set of genes, two approaches are available for identifying putative patterns of differential expression and regulation. A moment-based, biology-agnostic one uses a simple *t*-test to identify differences in the means of spliced counts in ``SearchData`` objects ``sd1`` and ``sd2``:

.. code-block:: python

 gf = compute_diffexp(sd1,sd2)
 
where ``gf`` is boolean vector that reports ``True`` if the gene is identified as DE. However, this approach cannot identify differences if biological parameters change in a correlated way and the mean stays the same. We introduce a more mechanistic approach for the identification of differential expression suggested by parameter variation, based on two ``SearchResults`` objects ``sr1`` and ``sr2``:

.. code-block:: python

 gf = compute_diffreg(sr1,sr2)
 
where ``gf`` is a two-dimensional boolean array that reports ``True`` if a particular *parameter* is identified as DE. After using these arrays to find a subpopulation of interest -- e.g., genes that do not exhibit variation in the spliced mean, but do exhibit modulation in the burst size -- it is possible to plug the gene filter ``genes_to_plot`` back in to inspect the raw data and fits:

.. code-block:: python

 gf = compare_gene_distributions(sr_arr,sd_arr,genes_to_plot=genes_to_plot)

