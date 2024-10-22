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
