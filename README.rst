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

Manual and examples
===================

Once installed, you can first check the version and input parameters with ``mquad -h`` 

MQuad recognizes 3 types of input: a cellSNP output folder (containing .vcf and AD/DP sparse matrices), AD/DP sparse matrix files (.mtx), or only the vcf file (cellSNP.cells.vcf.gz). Basic usage is as shown:

.. code-block:: bash

  mquad --vcfData $VCF -o $OUT_DIR -p 20