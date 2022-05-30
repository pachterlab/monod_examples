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
To generate an intronic index for `kb`, obtain a reference genome and run `kb ref`:

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
 
To then generate the count matrices in `loom` format, use `kb count`:

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

This generates a `loom` file in `OUTDIR/counts_filtered/`.

Pre-processing 
----------------

To define a "batch," or a set of inference runs, run :py:func:`preprocess.construct_batch`:

.. code-block:: python

 dir_string,dataset_strings = construct_batch(raw_filepaths, transcriptome_filepath, dataset_names, \
                                              attribute_names, batch_location, meta, batch_id, n_genes)

To import the files, specify the raw data (`loom`, `mtx`, or `adata`) filepaths in `raw_filepaths`, a gene length annotation filepath `raw_filepaths`, and various batch and dataset metadata. To specify the number of genes to analyze, set `n_genes`. 



//o retrieve a list of random ingredients,
//you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

