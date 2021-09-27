#Wrapper function for detecting useful mitochondrial variants

#import stuff
import os
from os import path
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import seaborn as sns
from scipy.io import mmread
from scipy.io import mmwrite
from scipy import sparse
from scipy.stats import betabinom, bernoulli, binom
import bbmix
from bbmix.models import MixtureBinomialSparseBatch
from kneed import KneeLocator
from collections import Counter
import vireoSNP
from .mquad_utils import findKnee

class MquadSparseMixBin():
    def __init__(self, AD, DP, variant_names=None, dataset_name=None):
        #initiate object with AD/DP sparse matrices
        #check if AD and DP have same length first
        self.ad = AD
        self.dp = DP

        if AD.shape[0] != DP.shape[0]:
            print('AD and DP length do not match!')
        else:
            print(str(AD.shape[0]) + ' variants detected...')

        if variant_names is not None:
            #sanity check for length of variant names
            if len(variant_names) != self.ad.shape[0]:
                print('No. of variant names does not match length of AD!')
            else:
                self.variants = variant_names
                print("Variant names detected...") 
            
        else:
            self.variants = None

        if dataset_name is not None:
            self.dataset = dataset_name
    
    def _check_outdir_exist(self, out_dir):
        if path.exists(out_dir) is not True:
            try:
                os.mkdir(out_dir)
                return True
            except:
                print("Can't make directory, do you have permission?")
                return False
        else:
            print('Out directory already exists, overwriting content inside...')
            return True
    
    def _addVariantNames(self, valid_rows):
        var_col = 'variant_name'
        variants = np.array(self.variants)
        df = self.df
        df_zeros = pd.DataFrame(np.zeros((variants.shape[0] - df.shape[0], df.shape[1])), columns=df.columns)
        df[var_col] = variants[valid_rows]
        df_zeros[var_col] = list(set(variants) - set(df[var_col]))
        self.df = pd.concat([df, df_zeros], axis=0, ignore_index=True)
        
    def _batchify(self, batch_size, x, y, valid_row_sizes):
        n_samples = valid_row_sizes.shape[0]
        start, item_start = 0, 0
        while start < n_samples:
            end = min(batch_size + start, n_samples)
            _row_sizes = valid_row_sizes[start:end]
            
            item_end = item_start + np.sum(_row_sizes)
            _x, _y = x[item_start:item_end], y[item_start:item_end]
            _ad = np.squeeze(np.asarray(self.ad[_x, _y]))
            _dp = np.squeeze(np.asarray(self.dp[_x, _y]))
            
            start = end
            item_start = item_end
            yield (_ad, _dp, _row_sizes)
        assert(item_start == x.shape[0])
            
    def validateSNP(self, SNP_names):
        ##input a list of SNP_names
        ##if the SNP has wrong reference, flag it

        ##read in ref vcf
        cell_vcf = vireoSNP.load_VCF("/home/aaronkwc/data/reproduced_results/cellsnp_ref/cellSNP.cells.vcf.gz", biallelic_only=True)
        variants=pd.Series(cell_vcf['variants'])
        df = pd.DataFrame(variants, columns=['variants'])
        df[['chr', 'pos', 'ref','alt']] = df['variants'].str.split('_', expand=True)
        #print(df.head())

        boo = []

        for name in SNP_names:
            pos, ref = name.split('_')[1], name.split('_')[2]
            #given a pos, find the ref and check if its correct
            #print(df.ref[df.pos == pos].values)
            if ((df.ref[df.pos == pos].values[0]) == ref) is True:
                boo.append(True)
            else:
                boo.append(False)

        return boo

    def fit_deltaBIC(self, out_dir, minDP=10, minAD=1, export_csv=True, nproc=30, batch_size=128):
        #here we fit and choose model based on deltaBIC
        n_variants = self.dp.shape[0]
        
        # adjust batch size to avoid unused processes
        adj_bs = min(n_variants // nproc, batch_size)
        print('CPUs used: {}, batch size: {} {}'.format(nproc, 
                                                        adj_bs, 
                                                        "" if adj_bs == batch_size else "(adjusted to avoid idle processes)"))
        print("Fitting in sparse mode...")
        t0=time.time()

        print("Initializing fit(mode: deltaBIC) on " + str(self.ad.shape[0]) + " variants...")
        
        dp_row, dp_col = np.nonzero(self.dp >= minDP)
        ad_row, _ = np.nonzero(self.ad >= minAD)
        
        # only variant with at leat one valid dp and ad records are included
        valid_rows = np.intersect1d(dp_row, ad_row)
        # filter invalid variants
        x, y = zip(*[(r, c) for r, c in zip(dp_row, dp_col) if r in valid_rows])
        
        # split batch
        x, y = np.array(x), np.array(y)
        valid_row_sizes = Counter(x)
        valid_row_sizes = np.array([valid_row_sizes[r_idx] for r_idx in valid_rows])
        
        assert(np.sum(valid_row_sizes) == x.shape[0])
        
        with mp.Pool(processes=nproc) as pool:
            results = pool.starmap_async(fit_batch, self._batchify(batch_size, x, y, valid_row_sizes)).get()
        self.df = pd.concat([pd.DataFrame(res) for res in results], axis=0, ignore_index=False)
        
        t1 = time.time()
        print("deltaBIC was calculated for " + str(self.ad.shape[0]) + " variants and took:%.2f minutes" %((t1-t0)/60))
        
#         self.df = pd.DataFrame()
#         for col, res in results.items():
#             self.df[col] = res.tolist()

        if self.variants is not None:
            self._addVariantNames(valid_rows)

        #sort df but keep index
        self.sorted_df = self.df.sort_values(by=['deltaBIC'], ascending=False)

        if export_csv is True:
            if self._check_outdir_exist(out_dir) is True:
                self.sorted_df.to_csv(out_dir + '/BIC_params.csv', index=False)
            else:
                self.sorted_df.to_csv('BIC_params.csv', index=False)

        self.df.to_csv(out_dir + '/debug_unsorted_BIC_params.csv', index=False)
        #return df of all metrics
        return self.df

    def selectInformativeVariants(self, min_cells=2, export_heatmap=True, export_mtx=True, out_dir=None):
        #takes self.df, return best_ad and best_dp as array

        if self.df is None:
            print('Fitted model not found! Have you run fit_deltaBIC/fit_logLik yet?')
        else:
            if out_dir is not None:
                if path.exists(out_dir) is not True:
                    try:
                        os.mkdir(out_dir)
                    except:
                        print("Can't make directory, do you have permission?")
            else:
                print('Out directory already exists, overwriting content inside...')

            x, y, knee, cutoff = findKnee(self.df.deltaBIC)
            
            plt.plot(x, y)
            plt.axvline(x=knee, color="black", linestyle='--',label="cutoff")
            plt.legend()
            plt.ylabel("\u0394BIC")
            plt.xlabel("Cumulative probability")
            plt.savefig(out_dir + '/deltaBIC_cdf.pdf')

            #make a PASS/FAIL column in self.df for easier subsetting
            print(cutoff)
            #self.sorted_df['VALID'] = self.validateSNP(self.sorted_df.variant_name)
            self.sorted_df['PASS_KP'] = self.sorted_df.deltaBIC.apply(lambda x: True if x >= cutoff else False)
            self.sorted_df['PASS_MINCELLS'] = self.sorted_df.num_cells_minor_cpt.apply(lambda x: True if x >= min_cells else False)

            self.final_df = self.sorted_df[(self.sorted_df.PASS_KP == True) & (self.sorted_df.PASS_MINCELLS == True)]
            #print(self.final_df.head())
            
            #will deprecate in later versions
            #self.final_df = self.sorted_df[0:int(len(y) * (1 - knee))]
            #self.final_df = self.final_df[self.sorted_df.num_cells_minor_cpt >= min_cells]

            print('Number of variants passing threshold: '  + str(len(self.final_df['variant_name'])))

            if len(self.final_df['variant_name']) != 0:
                passed_variants = self.final_df['variant_name']
                idx = [self.variants.index(i) for i in passed_variants]
        
                best_ad = self.ad[idx]
                best_dp = self.dp[idx]
            else:
                print("No informative variants detected! If you are using 10x data, try setting --minDP to a smaller number.")


        self.sorted_df.to_csv(out_dir + '/BIC_params.csv', index=False)
        #fname = by + '_' + str(threshold) + '_'

        if self.variants is not None:
            #best_vars = np.array(self.variants)[idx]
            renamed_vars = []
            for var in passed_variants:
                renamed_vars.append((var.split('_')[1] + var.split('_')[2] + '>' + var.split('_')[3]))

            with open(out_dir + '/passed_variant_names.txt', "w+") as var_file:
                var_file.write('\n'.join(str(var) for var in renamed_vars))
                
        if export_heatmap is True:
            af = best_ad/best_dp
            #print(af.shape)
            #af = af.fillna(0)
            fig, ax = plt.subplots(figsize=(15,10))
            plt.title("Allele frequency of top variants")
            plt.style.use('seaborn-dark')
            if self.variants is not None:
                sns.heatmap(af, cmap='Greens', yticklabels=renamed_vars)
            else:
                sns.heatmap(af, cmap='Greens')
            plt.savefig(out_dir + '/top variants heatmap.pdf')

        #export ad dp mtx out for vireo
        if export_mtx is True:
            mmwrite(out_dir + '/passed_ad.mtx', sparse.csr_matrix(best_ad))
            mmwrite(out_dir + '/passed_dp.mtx', sparse.csr_matrix(best_dp))

        return best_ad, best_dp

    def readParams(self, file):
        self.df = pd.read_csv(file)
        self.sorted_df = self.df.sort_values(by=['deltaBIC'], ascending=False)

        return self.df, self.sorted_df
    
def sparseMixBinFit(valid_ad, valid_dp, valid_row_sizes, fix_seed=None):
    #input ad dp arrays, output params, BICs, delta BIC        
    if fix_seed is not None:
        np.random.seed(fix_seed)

    model1 = MixtureBinomialSparseBatch(n_components = 1, tor=1e-20)
    params1 = model1.fit((valid_ad, valid_dp), valid_row_sizes, max_iters=500, early_stop=True)

    model2 = MixtureBinomialSparseBatch(n_components = 2,tor=1e-20)
    params2 = model2.fit((valid_ad, valid_dp), valid_row_sizes, max_iters=500, early_stop=True)

    delta_BIC = model1.model_scores["BIC"] - model2.model_scores["BIC"]

    p = params2[:, [0,1]]
    pi = params2[:, [2,3]]
    fraction_b_allele = np.min(p, axis=1) * np.array([pi[ith, idx] for ith, idx in enumerate(np.argmin(p, axis=1))])

    new_mutation = np.zeros(p.shape[0], dtype=bool)
    as_mutation = np.zeros(p.shape[0], dtype=bool)

    # new_mutation
    new_mut_sel = (np.max(pi, axis=1) < 0.95) & (np.min(p, axis=1) < 0.05) & (np.max(p, axis=1) > 0.1)
    new_mutation[new_mut_sel] = True

    # as_mutation
    as_mut_sel = (np.min(p, axis=1) > 0.1) & (np.min(pi, axis=1) > 0.15) & (~new_mut_sel)
    as_mutation[as_mut_sel] = True

    minor_cpt_n = np.min(pi, axis=1) * valid_row_sizes

    results = {
        "num_cells": valid_row_sizes,
        'deltaBIC': delta_BIC,
        'params1': params1.tolist(),
        'params2': params2.tolist(),
        'model1BIC': model1.model_scores["BIC"],
        'model2BIC': model2.model_scores["BIC"],
        'new_mutation': new_mutation, 
        'as_mutation': as_mutation, 
        'fraction_b_allele': fraction_b_allele, 
        'num_cells_minor_cpt': minor_cpt_n, 
    }
    return results    

def fit_batch(valid_ad, valid_dp, valid_row_sizes):
    basic_stats = basicStats(valid_ad, valid_dp, valid_row_sizes)
    results = sparseMixBinFit(valid_ad, valid_dp, valid_row_sizes)
    results.update(basic_stats)
    return results
    
def basicStats(valid_ad, valid_dp, valid_row_sizes):
    #basic staistics
    #Total DP across all cells
    left, batch_size = 0, len(valid_row_sizes)
    stats = ['total_DP', 'median_DP', 'total_AD', 'median_AD', 'num_cells_nonzero_AD']
    batch_res = {name:np.empty(batch_size) for name in stats}

    for ith, smp_sz in enumerate(valid_row_sizes):
        right = left + smp_sz
        _d = valid_dp[left:right]
        _a = valid_ad[left:right]
        batch_res[stats[0]][ith] = np.sum(_d)
        batch_res[stats[1]][ith] = np.median(_d)
        batch_res[stats[2]][ith] = np.sum(_a)
        batch_res[stats[3]][ith] = np.median(_a)
        batch_res[stats[4]][ith] = np.count_nonzero(_a)

        left = right

    return batch_res

    total_DP = np.sum(_d, axis=1)
    #Median DP across all cells
    median_DP = np.median(_d, axis=1)
    #Total AD across all cells
    total_AD = np.sum(_a, axis=1)
    #Median AD across all cells
    median_AD = np.median(_a, axis=1)
    #How many cells have this variant?
    non_zero = np.count_nonzero(_a, axis=1)
    return {'total_DP' :total_DP,
            'median_DP':median_DP,
            'total_AD' :total_AD,
            'median_AD':median_AD,
            'num_cells_nonzero_AD':non_zero
           }    

if __name__ == '__main__':
    import vireoSNP
    from vireoSNP.utils.io_utils import read_sparse_GeneINFO
    from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb
    
    cell_vcf = vireoSNP.load_VCF("example/example.vcf.gz", biallelic_only=True)
    cell_dat = vireoSNP.vcf.read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])
    mdphd = MquadSparseMixBin(AD = cell_dat['AD'], DP = cell_dat['DP'], variant_names= cell_vcf['variants'])

    df = mdphd.fit_deltaBIC(out_dir='test') 
    mdphd.selectInformativeVariants(out_dir = 'test')
