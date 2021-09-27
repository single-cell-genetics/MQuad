# CLI for mitoMut

import os
import sys
import time
import numpy as np
from scipy.io import mmread
from optparse import OptionParser, OptionGroup

from .version import __version__
from vireoSNP.utils.io_utils import read_cellSNP, read_vartrix, read_sparse_GeneINFO
from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb

from .mquad import Mquad
from .mquad_batch_mixbin import MquadSparseMixBin

START_TIME = time.time()

def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--cellData", "-c", dest="cell_data", default=None,
        help=("The cellSNP folder with AD and DP sparse matrices."))
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
        help=("Dirtectory for output files [default: $cellData/mitoMut]"))
    
    group0 = OptionGroup(parser, "Alternative input formats")
    group0.add_option("--mtxData", "-m", dest="mtx_data", default=None,
        help=("The two mtx files for AD and DP matrices, comma separated"))
    parser.add_option("--vcfData", dest="vcf_data", default=None,
        help=("The cell genotype file in VCF format"))
    group0.add_option("--BICparams", "--b", dest="BIC_params", default=None,
        help=("Existing unsorted_debug_BIC_params.csv"))
    group0.add_option("--tenx", "--t", dest="cutoff", default=None,
        help=("User-defined deltaBIC cutoff mainly for low-depth data"))

    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option("--randSeed", type="int", dest="rand_seed", default=None,
        help="Seed for random initialization [default: %default]")
    group1.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")
    group1.add_option("--minDP", type="int", dest="minDP", default=10, 
        help="Minimum DP to include for modelling [default: 10]")
    group1.add_option("--minCell", type='int', dest="minCell", default=2,
        help=("Minimum no. of cells in minor component [default: 2]"))
    group1.add_option("--batchFit", type='int', dest="batch_fit", default=1,
        help=("1 if fit MixBin model using batch mode, 0 else [default: 1]"))
    group1.add_option("--batchSize", type='int', dest="batch_size", default=128,
        help=("Number of variants in one batch, cooperate with --nproc for speeding up [default: 128]"))
    
    
    parser.add_option_group(group0)
    parser.add_option_group(group1)
    (options, args) = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        print("Welcome to MQuad v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)

    ## out directory
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $cellFilePath/mquad")
        out_dir = os.path.dirname(os.path.abspath(options.cell_file)) + "/mquad"
    elif os.path.dirname(options.out_dir) == "":
        out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ## input data (VCF.gz or a folder with sparse matrices)
    if (options.cell_data is None and options.mtx_data is None and 
        options.vcf_data is None):
        print("Error: need cell data in vcf file, or cellSNP output folder, or "
              "matrix files: AD.mtx,DP.mtx")
        sys.exit(1)
    elif options.mtx_data is not None:
        print("[MQuad] Loading matrix files ...")
        matrix_files = options.mtx_data.split(",")
        if len(matrix_files) != 2:
            print("Error: mtxData requires 2 comma separated files")
            sys.exit(1)
        cell_dat = {}
        cell_dat['AD'] = mmread(matrix_files[0])
        cell_dat['DP'] = mmread(matrix_files[1])
        cell_dat['variants'] = ['SNP%d' %(x + 1) for x in 
                                range(cell_dat['AD'].shape[0])]
        
    elif options.vcf_data is not None:
        print("[MQuad] Loading cell VCF file ...")
        cell_vcf = load_VCF(options.vcf_data, biallelic_only=True)
        cell_dat = read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
        for _key in ['samples', 'variants', 'FixedINFO', 'contigs', 'comments']:
            cell_dat[_key] = cell_vcf[_key]
        
    else:
        #os.path.isdir(os.path.abspath(options.cell_data)):
        print("[MQuad] Loading cell folder ...")
        cell_dat = read_cellSNP(options.cell_data)
        
    ## More options
    nproc = options.nproc
    minDP = options.minDP
    batch_size = options.batch_size
    cutoff = options.cutoff
    minCell = options.minCell
    
    ## Main functions
    if options.BIC_params is not None:
        mdphd = Mquad(AD = cell_dat['AD'], DP = cell_dat['DP'], 
                        variant_names = cell_dat['variants'])
        print("[MQuad] Using existing BIC params to filter variants only...")
        best_ad, best_dp = mdphd.selectInformativeVariants(min_cells = minCell, out_dir = out_dir, existing_df=options.BIC_params, tenx_cutoff=cutoff)
    else:
        if options.batch_fit == 0:
            mdphd = Mquad(AD = cell_dat['AD'], DP = cell_dat['DP'], 
                            variant_names = cell_dat['variants'])
            df = mdphd.fit_deltaBIC(out_dir = out_dir, nproc = nproc, minDP = minDP, beta_mode = False)
            best_ad, best_dp = mdphd.selectInformativeVariants(min_cells = minCell, out_dir = out_dir, tenx_cutoff=cutoff)
        else:
            #use sparse mode for faster performance
            #default to sparse mode in v1.6.0
            mdphd = MquadSparseMixBin(
                AD=cell_dat['AD'], 
                DP=cell_dat['DP'], 
                variant_names=cell_dat['variants']
            )
            df = mdphd.fit_deltaBIC(
                out_dir=out_dir,
                minDP=minDP,
                nproc=nproc,
                batch_size=batch_size
            )
            best_ad, best_dp = mdphd.selectInformativeVariants(min_cells = minCell, out_dir = out_dir)
    
    run_time = time.time() - START_TIME
    print("[MQuad] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                  run_time % 60))
    print()
    
        
if __name__ == "__main__":
    main()

