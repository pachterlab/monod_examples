# import packages
import os,sys
sys.path.insert(0, '../monod/src/')

import numpy as np
import pickle
import time
import anndata

# import Monod
import monod
from monod import preprocess, extract_data, cme_toolbox, inference, analysis


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

gradient_params = {'max_iterations':20,'init_pattern':'moments','num_restarts':1}

phys_lb = [-1.0, -1.8, -1.8 ]
phys_ub = [4.2, 2.5, 3.5]
samp_lb = [-9, -4] #-7.5, -2
samp_ub = [-4, 1.5] #-5.5, 0
gridsize = [10,11]


timing_dict = { samp : [] for samp in samples }

loom_filepaths = [f'./data/lu_radiation/{sample}/counts_unfiltered/adata_T.loom' for sample in samples] 
n_datasets = len(loom_filepaths)
dir_string,dataset_strings = monod.preprocess.construct_batch(loom_filepaths,                                        transcriptome_filepath,dataset_names,attribute_names=attribute_names,                                             batch_location='./fits',meta=f'lu_T',batch_id=1,creator='mc',                                     exp_filter_threshold=None,cf=None,viz=False,n_genes=3000,filt_param=filt_param,datestring='241031')

result_strings = []

        # run inference
for i in range(n_datasets):
    try:
        t1 = time.time()
        fitmodel = monod.cme_toolbox.CMEModel('Bursty','Poisson')
        inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,                dataset_strings[i],fitmodel,use_lengths = True,
            gradient_params = gradient_params)
        search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
            dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=None)
        full_result_string = inference_parameters.fit_all_grid_points(20,search_data)
        result_strings.append(full_result_string)
        t2 = time.time()
        timing_dict[samples[i]] = t2-t1
    except:
        continue
        
# save timing dictionary
with open('results/lu_T_timing.pickle', 'wb') as handle:
    pickle.dump(timing_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)