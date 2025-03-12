from memory_profiler import profile

# import packages
import os,sys
sys.path.insert(0, '../monod/src/')

import numpy as np
import loompy as lp
import pandas as pd

# timing
import time

# import Monod
import monod
from monod import preprocess, extract_data, cme_toolbox, inference, analysis


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-n_genes","--n_genes",type=int,default=1000)
parser.add_argument("-sample","--sample",type=int,default=0)
args = parser.parse_args()
n_genes = args.n_genes
s = args.sample


num_cores = 20

gradient_params = {'max_iterations':2,'init_pattern':'moments','num_restarts':1}

@profile        
def run(n_genes=n_genes):
    

    samples = ['S0',
        'S22',
        'S23',
        'S26',
    #     'S1',
        'S2',
    #     'S3',
        'S35',
        'S10',
        'S8',
        'S9',
        'S14',
        'S15',
        'S17',
    #     'S37',
        'S38',
        'S40']
    
    samples = samples[s:s+1]

    dataset_meta = [samp + '_T' for samp in samples]
    dataset_names = dataset_meta 

    transcriptome_filepath = './transcriptome/gg_200525_genome_polyA_cum_3.txt'

    spliced_layer = 'spliced'
    unspliced_layer = 'unspliced'
    gene_attr = 'GENE_NAME'
    cell_attr = 'barcode'

    attribute_names=[(unspliced_layer,spliced_layer),gene_attr,cell_attr]

    filt_param={
            "min_U_mean": 0.01,
            "min_S_mean": 0.01,
            "max_U_max": 1000,
            "max_S_max": 1000,
            "min_U_max": 3,
            "min_S_max": 3,
        }

    phys_lb = [-1.0, -1.8, -1.8 ]
    phys_ub = [4.2, 2.5, 3.5]
    samp_lb = [-9, -4] #-7.5, -2
    samp_ub = [-4, 1.5] #-5.5, 0
    gridsize = [1,1]
    
    i = 0

    loom_filepaths = [f'./data/lu_radiation/{sample}/counts_unfiltered/adata_T.loom' for sample in samples] 
    n_datasets = len(loom_filepaths)
    dir_string,dataset_strings = monod.preprocess.construct_batch(loom_filepaths,                                        transcriptome_filepath,dataset_names,attribute_names=attribute_names,                                             batch_location='./fits',meta=f'lu_T',batch_id=1,creator='mc',                                     exp_filter_threshold=None,cf=None,viz=False,n_genes=n_genes,filt_param=filt_param,datestring='241212')

    result_strings = []

    # run inference

    fitmodel = monod.cme_toolbox.CMEModel('Bursty','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = gradient_params)
    search_data1a = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
        dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=None)
    full_result_string = inference_parameters.fit_all_grid_points(20,search_data1a)
    result_strings.append(full_result_string)

        
if __name__ == "__main__" :
    run()


