=====
MQuad
=====

MQuad: Mixture Model for Mitochondrial Mutation detection in single-cell omics data

MQuad is a tool that detects mitochondrial mutations that are informative for clonal substructure inference. It uses a binomial mixture model to assess the heteroplasmy of mtDNA variants among background noise.

A recommended pipeline to generate the neccessary files:

1. use `cellSNP <https://github.com/single-cell-genetics/cellSNP>`_ or `cellsnp-lite <https://github.com/single-cell-genetics/cellsnp-lite>`_ (a faster version of cellSNP, still at testing stage so might be unstable) to pileup mtDNA variants from raw .bam file(s)

2. use MQuad to differentiate informative mtDNA variants from noisy backbground

3. use `vireoSNP <https://github.com/single-cell-genetics/vireo>`_ to assign cells to clones based on mtDNA variant profile

Different upstream/downstream packages can also be used if the neccesary file formats are available.

Installation
============

MQuad is available through Python Package Index. To install, type the following command line and add ``-U`` for updates:

.. code-block:: bash

  pip install -U mquad

Alternatively, you can install from this GitHub repository for latest (often development) version by the following command line:

.. code-block:: bash

  pip install -U git+https://github.com/single-cell-genetics/MQuad

Manual
======

Once installed, you can first check the version and input parameters with ``mquad -h`` 

MQuad recognizes 3 types of input:

1. cellSNP output folder with AD and DP sparse matrices (.mtx)

.. code-block:: bash

  mquad -c $INPUT_DIR -o $OUT_DIR -p 20

2. .vcf only

.. code-block:: bash

  mquad --vcfData $VCF -o $OUT_DIR -p 20

3. AD and DP sparse matrices (.mtx), comma separated

.. code-block:: bash

  mquad -m cellSNP.tag.AD.mtx, cellSNP.tag.DP.mtx -o $OUT_DIR -p 20
  
For droplet-based sequencing data, eg. 10X Chromium CNV, scATAC..etc, it is recommended to add ``--minDP 5`` or a smaller value to prevent errors during fitting. The default value is 10, which is suitable for Smart-seq2 data but might be too stringent for low sequencing depth data.

The output files will be explained below in the 'Example' section.

Example
=======

MQuad comes with an example dataset for you to test things out. The mtDNA mutations of this dataset are extracted from `Ludwig et al, Cell, 2019 <https://doi.org/10.1016/j.cell.2019.01.022>`_. It contains 500 background variants, along with 9 variants used in Supp Fig. 2F (and main Fig. 2F). There is also 1 additional variant that is informative but not mentioned in the paper. In total, there are 510 variants in the example dataset.

Run the following command line:

.. code-block:: bash

  mquad --vcfData example/example.vcf.gz -o example_test -p 5
  
The output files should include:

* passed_ad.mtx, passed_dp.mtx: Sparse matrix files of the AD/DP of qualified variants for downstream clonal analysis
* top variants heatmap.pdf: Heatmap of the allele frequency of qualified variants

.. image:: images/top_var_heatmap.png
    :width: 100px
    :align: center
    :height: 50px
    
* deltaBIC_cdf.pdf: A cdf plot of deltaBIC distribution of all variants, including the cutoff determined by MQuad

.. image:: images/cdf.png
    :width: 100px
    :align: center
    :height: 50px
    
* BIC_params.csv: A spreadsheet containing detailed parameters/statistics of all variants, sorted from highest deltaBIC to lowest
* debug_unsorted_BIC_params.csv: Same spreadsheet as BIC_params.csv but unsorted, for developers' debugging purpose, will probably be removed on later versions of MQuad
