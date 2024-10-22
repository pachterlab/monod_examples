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
