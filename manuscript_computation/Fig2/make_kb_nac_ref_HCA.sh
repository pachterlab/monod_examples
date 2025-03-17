mkdir /home/tchari/metadata/refdata-gex-GRCh38-2020-A/

mkdir /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28

kb ref -i /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/index.idx \
-g /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/t2g.txt \
-f1 /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/cdna.fa \
-f2 /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/nascent.fa \
-c1 /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/cdna_t2c.txt \
-c2 /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/nascent_t2c.txt \
--workflow nac \
--overwrite \
/home/ggorin/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
/home/ggorin/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf
