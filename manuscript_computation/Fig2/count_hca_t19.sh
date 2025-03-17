kb count --verbose -i /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/index.idx \
-g /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/t2g.txt -x 10xv3 \
-o ./T19 -t 10 -m 30G \
-c1 /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/cdna_t2c.txt \
-c2 /home/tchari/metadata/refdata-gex-GRCh38-2020-A/kallisto_nac_28/nascent_t2c.txt \
--workflow nac --sum cell --overwrite --loom \
/home/tchari/counts/monod_rev_HCA/fastq/SRR19039313_2.fastq /home/tchari/counts/monod_rev_HCA/fastq/SRR19039313_3.fastq \
/home/tchari/counts/monod_rev_HCA/fastq/SRR19039314_2.fastq /home/tchari/counts/monod_rev_HCA/fastq/SRR19039314_3.fastq
