Usage
=====

.. _installation:

Installation
------------

To use *Monod*, install it from `pip`:

.. code-block:: console

 pip install monod

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
 
To then generate the count matrices in ``loom`` format, use ``kb count``:

.. code-block:: console

   kb count --verbose 
   -i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx 
   -g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt 
   -x 10xv3 
   -o OUTDIR/ 
   -t 30 -m 30G 
   -c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt 
   -c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt 
   --workflow lamanno --filter bustools --overwrite --loom 
   DATASET_FASTQ_LOCATIONS

This generates a ``loom`` file in ``OUTDIR/counts_filtered/``.

Pre-processing 
----------------

To define a "batch," or a set of inference runs, run :py:func:`preprocess.construct_batch`:

.. code-block:: python

 dir_string,dataset_strings = construct_batch(raw_filepaths, transcriptome_filepath, dataset_names, \
                                              attribute_names, batch_location, meta, batch_id, n_genes)

To import the files, specify the raw data (``loom``, ``mtx``, or ``adata``) filepaths in ``raw_filepaths``, a gene length annotation filepath ``transcriptome_filepath``, and various batch and dataset metadata. To specify the number of genes to analyze, set ``n_genes``. 

This will create a batch directory, dataset-specific subdirectories, the file ``gene_set.csv`` with the list of genes that meet filtering thresholds for all datasets, and the file ``genes.csv`` with the list of genes selected for further analysis. 

Model, data, and parameter definition 
----------------

To define a CME model of transcription and sequencing, initialize an instance of :py:class:`cme_toolbox.CMEModel`:

.. code-block:: python

 CMEModel(biological_model,sequencing_model)

where ``biological_model = {'Bursty','Constitutive','Extrinsic','CIR'}`` represents the transcriptional process and ``sequencing_model = {'None','Bernoulli','Poisson'}`` represents the dynamics of the sampling process.

To define the search parameters, initialize an instance of :py:class:`inference.InferenceParameters`:

.. code-block:: python
 inference_parameters = InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                     dataset_string,fitmodel,use_lengths)

where ``phys_lb`` and ``phys_ub`` are bounds on the transcriptional process model parameters, ``samp_lb`` and ``samp_ub`` are bounds on the sampling process model parameters, ``gridsize`` defines the grid for the sampling parameter scan, and ``use_lengths`` determines whether the unspliced mRNA capture rate depends on the gene length (to model priming at ubiquitous internal polyA sites).

To create a ``SearchData`` object to input into the inference process, run :py:func:`extract_data.extract_data`:

.. code-block:: python
 extract_data(loom_filepath, transcriptome_filepath, dataset_name,
                                dataset_string, dir_string)

Running the inference pipeline 
----------------

To run the pipeline, simply call the following parallelized code:

.. code-block:: python
 inference_parameters.fit_all_grid_points(n_cores,search_data)

