Usage
=====

Pre-processing 
----------------

In the anndata Monod version, cell/gene filtering can be done using scanpy, and the resulting anndata object can be directly input into the inference procedure. For example:

.. code-block:: python

 gaba_adata = anndata.read_h5ad('./gaba_example.h5ad')
 scanpy.pp.filter_cells(gaba_adata, min_counts=threshold)

..
 add more here?


Defining a Model
----------------------

To define a CME model of transcription and sequencing, initialize an instance of :py:class:`cme_toolbox.CMEModel`:

.. code-block:: python

 fitmodel = cme_toolbox.CMEModel('ProteinBursty','Poisson')

where ``biological_model = {'Bursty','Constitutive','Extrinsic','CIR'}`` represents the transcriptional process and ``sequencing_model = {'None','Bernoulli','Poisson'}`` represents the dynamics of the sampling process.


Inference
----------------

The entire inference procedure can be carried out in a single function, :py:func:`inference.perform_inference`, using anndata or an h5ad filename as input.

.. code-block:: python

 fitted_adata = inference.perform_inference(input_adata, fitmodel, dataset_string='your_output_dirname')

This call will first filter the given anndata for suitable genes. Then it will iterate over a grid of sampling parameter values, and fit the optimum biological parameters for each gene for the given CME model.

The output is an anndata object including moments and fitted parameters for each gene, at the sampling parameter optimum. The results will be saved to the output directory specified in ``dataset_string``. 

Length-dependent noise models can be implemented if a transcriptome length annotation is provided, using ``transcriptome_filepath='transcriptome_filepath'``.

Fitting Parameters by Cluster
-------------------

Alternatively one can define the search parameters and cluster the data. This will run the `meK-Means <https://github.com/pachterlab/CGP_2023/>`_ clustering algorithm (see the paper `here <https://www.biorxiv.org/content/10.1101/2023.09.17.558131v2>`_). In this case, the user should give a value for ``mek_means_params = (k, epochs)``, where ``k`` is the user-defined number of clusters to learn and ``epochs`` is the numbers of rounds to learn the data clusters. All other parameters remain the same. 

Post-processing and QC
----------------

To visualize the results of the fitting procedure, you can run: 

.. code-block:: python

 analysis.run_qc(fitted_adata)

This will plot a color map of the fitted distribution for each gene against its observed count distribution.

It will also visualize the stability of the fits under subsampling, and the length dependence of the fitted parameters.



TODO: Do we want to add noise decomposition?

Differential parameter value identification
----------------

Given two fitted anndata objects, we can analyze the differential parameters between the two datasets, as long as some genes overlap.

We can run: 

.. code-block:: python

 DE_genes, DE_filter, offs, residuals = analysis.DE_parameters(fitted_adata_1, fitted_adata_2)

This will output a list of genes with signficantly different parameters between datasets, along with their offsets and residuals (TODO: explain).

It will also modify the anndata objects, adding columns for the fold-changes in parameter values between genes in the two datasets.

If we have fitted using meK-Means, we can perform differential parameter analysis between clusters in the same way, using just one anndata object:

.. code-block:: python

 DE_genes, DE_filter, offs, residuals = analysis.DE_parameters(fitted_adata_mek)
 
