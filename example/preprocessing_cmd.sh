#This file describes pre-processing cellsnp-lite commands before using MQuad, assuming you have everything installed correctly

### SMART-SEQ2 ####
#for sequencing data with individual .bam files for each cell + no UMItags/barcodes

#run cellsnp-lite mode2a on bam list
#change --chrom= to whatever reference genome you aligned to - in this case we use hg19
ls *.bam > bam.lst
cellsnp-lite -S bam.lst -i bam.lst -o cellsnp --UMItag None --genotype --gzip --chrom=chrM -p 10


### 10X/UMI-BASED ###
#for 10x and UMI-based sequencing data where there is only 1 big .bam + barcodes

#run cellsnp-lite on the .bam directly
#change --chrom= to whatever reference genome you aligned to, in most cases 10x data are aligned to GrCh38 so the chr name is MT
cellsnp-lite -s possorted_genome.bam -b barcodes.tsv -o cellsnp --chrom=MT --UMItag Auto --minMAF 0 --minCOUNT 0 --genotype --gzip -p 10


#The above steps should generate a cell x snp .vcf file (cellSNP.cells.vcf.gz), or AD/DP sparse matrices if you did not use the --genotype option

#Then you should be able to run MQuad on the vcf file 

#More tutorials and updates coming...
