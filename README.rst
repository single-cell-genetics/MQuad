=====
MQuad
=====

MQuad: Mixture Model for Mitochondrial Mutation detection in single-cell omics data

MQuad is a tool that detects mitochondrial mutations that are informative for 
clonal substructure inference. 

Installation
============

MQuad is available through PyPI_. To install, type the following command line and add ``-U`` for updates:

.. code-block:: bash

  pip install -U mquad

Alternatively, you can install from this GitHub repository for latest (often development) version by the following command line:

.. code-block:: bash

  pip install -U git+https://github.com/single-cell-genetics/MQuad

Manual
======

Once installed, you can first check the version and input parameters with ``mquad -h`` 

MQuad recognizes 3 types of input: a cellSNP output folder (containing .vcf and AD/DP sparse matrices), AD/DP sparse matrix files (.mtx), or only the vcf file (cellSNP.cells.vcf.gz). Basic usage is as shown:

.. code-block:: bash

  mquad --vcfData $VCF -o $OUT_DIR -p 20
  
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
    :width: 200px
    :align: center
    :height: 100px
    
* deltaBIC_cdf.pdf: A cdf plot of deltaBIC distribution of all variants, including the cutoff determined by MQuad

.. image:: images/cdf.png
    :width: 200px
    :align: center
    :height: 100px
    
* BIC_params.csv: A spreadsheet containing detailed parameters/statistics of all variants, sorted from highest deltaBIC to lowest
* debug_unsorted_BIC_params.csv: Same spreadsheet as BIC_params.csv but unsorted, for developers' debugging purpose, will probably be removed on later versions of MQuad
