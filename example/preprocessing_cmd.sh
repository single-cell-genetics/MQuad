#This file describes pre-processing cellsnp-lite commands before using MQuad, assuming you have everything installed correctly

### SMART-SEQ2 ####
#for sequencing data with individual .bam files for each cell + no UMItags/barcodes

#first use mode2 on a merged bulk .bam
ls *.bam > bam.lst
samtools merge -b bam.lst out.bam
samtools index out.bam

#generate .vcf file for all cells in bulk
cellsnp-lite -s out.bam -o mode2 --chrom=MT --UMItag None --minMAF 0 --minCOUNT 0 --genotype --gzip -p 10

#pileup for each cell with the previous bulk .vcf as input
cellsnp-lite -S bam.lst -I bam.lst -o mode3 -R mode2/cellSNP.cells.vcf.gz --UMItag None --genotype --gzip -p 10


### 10X/UMI-BASED ###
#for 10x and UMI-based sequencing data where there is only 1 big .bam + barcodes

#run mode2 on the .bam directly
cellsnp-lite -s possorted_genome.bam -b barcodes.tsv -o mode2 --chrom=MT --UMItag Auto --minMAF 0 --minCOUNT 0 --genotype --gzip -p 10


#The above steps should generate a cell x snp .vcf file (cellSNP.cells.vcf.gz), or AD/DP sparse matrices if you did not use the --genotype option

#More tutorials and updates coming...
