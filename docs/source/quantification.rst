Quantification
=====

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
