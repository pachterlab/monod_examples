Quantification
=====

To generate nascent and mature count matrices for analysis with Monod, use `kb-python`_, using the option ``--h5ad`` in ``kb count`` to generate anndata output, and ``--sum=cell`` to add mature and ambiguous counts into a single matrix.

.. _kb-python: https://www.nature.com/articles/s41596-024-01057-0

To generate an intronic index for ``kb``, obtain a reference genome and run ``kb ref``:

.. code-block:: console

 kb ref 
 --workflow=nac 
 -i index.idx 
 -g t2g.txt 
 -c1 cdna.txt 
 -c2 nascent.txt
 -f1 cdna.fasta 
 -f2 nascent.fasta 
 genome.fasta
 genome.gtf
 
To then generate the count matrices in ``h5ad`` format, use ``kb count``:

.. code-block:: console

 kb count 
 -x <tech> 
 --workflow=nac 
 -w onlist.txt
 -o output_dir 
 -i index.idx 
 -g t2g.txt
 -c1 cdna.txt 
 -c2 nascent.txt 
 --sum=cell 
 --h5ad 
 R1.fastq
 R2.fastq

This generates an ``h5ad`` file in ``OUTDIR/counts_filtered/``.
