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

 fitted_adata = inference.perform_inference(input_adata, fitmodel, dataset_string='your_output_dirname')

This call will first filter the given anndata for suitable genes. Then it will iterate over a grid of sampling parameter values, and fit the optimum biological parameters for each gene for the given CME model.

The output is an anndata object including moments and fitted parameters for each gene, at the sampling parameter optimum. The results will be saved to the output directory specified in ``dataset_string``. 

Length-dependent noise models can be implemented if a transcriptome length annotation is provided, using ``transcriptome_filepath='transcriptome_filepath'``.

Fitting Parameters by Cluster
-------------------

Alternatively one can define the search parameters and cluster the data. This will run the `meK-Means <https://github.com/pachterlab/CGP_2023/>`_ clustering algorithm (see the paper `here <https://www.biorxiv.org/content/10.1101/2023.09.17.558131v2>`_). In this case, the user should give a value for ``mek_means_params'' = (``k``, ``epochs``), where ``k`` is the user-defined number of clusters to learn and ``epochs`` is the numbers of rounds to learn the data clusters. All other parameters remain the same. 

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
Given a set of matched datasets, run with the same model over the same set of genes, two approaches are available for identifying putative patterns of differential expression and regulation. A moment-based, biology-agnostic one uses a simple *t*-test to identify differences in the means of spliced counts in ``SearchData`` objects ``sd1`` and ``sd2``:

.. code-block:: python

 gf = compute_diffexp(sd1,sd2)
 
where ``gf`` is boolean vector that reports ``True`` if the gene is identified as DE. However, this approach cannot identify differences if biological parameters change in a correlated way and the mean stays the same. We introduce a more mechanistic approach for the identification of differential expression suggested by parameter variation, based on two ``SearchResults`` objects ``sr1`` and ``sr2``:

.. code-block:: python

 gf = compute_diffreg(sr1,sr2)
 
where ``gf`` is a two-dimensional boolean array that reports ``True`` if a particular *parameter* is identified as DE. After using these arrays to find a subpopulation of interest -- e.g., genes that do not exhibit variation in the spliced mean, but do exhibit modulation in the burst size -- it is possible to plug the gene filter ``genes_to_plot`` back in to inspect the raw data and fits:

.. code-block:: python

 gf = compare_gene_distributions(sr_arr,sd_arr,genes_to_plot=genes_to_plot)
 
