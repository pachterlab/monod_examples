# import packages
import os,sys
sys.path.insert(0, '../monod/src/')

import numpy as np
import loompy as lp
import anndata
import pandas as pd

# colors
import matplotlib.pyplot as plt

# import Monod
import monod
from monod import preprocess, extract_data, cme_toolbox, inference, analysis


dataset_meta = ['mESC',]
dataset_names = dataset_meta 
raw_data_locations = dataset_meta

transcriptome_filepath = './transcriptome/gg_200525_genome_polyA_cum_3.txt'

spliced_layer = 'spliced'
unspliced_layer = 'unspliced'
gene_attr = 'GENE_NAME'
cell_attr = 'cell_bc'

attribute_names=[(unspliced_layer,spliced_layer),gene_attr,cell_attr]


data_path = './data/khateb/'

loom_filepaths = [data_path+x+'.loom' for x in raw_data_locations] 
n_datasets = len(loom_filepaths)


# # Preprocess
import logging, sys
logging.basicConfig(stream=sys.stdout)
log = logging.getLogger()
log.setLevel(logging.INFO)


#Filter again for acceptable U/S bounds

filt_param={
        "min_U_mean": 0.01,
        "min_S_mean": 0.01,
        "max_U_max": 1000,
        "max_S_max": 1000,
        "min_U_max": 1,
        "min_S_max": 1,
    }



dir_string,dataset_strings = monod.preprocess.construct_batch(loom_filepaths,transcriptome_filepath,dataset_names,attribute_names=attribute_names,batch_location='./fits',meta='khateb',batch_id=1,creator='mc',exp_filter_threshold=None,cf=None,viz=False,n_genes=2000,filt_param=filt_param)


# # Run inference

gradient_params = {'max_iterations':20,'init_pattern':'moments','num_restarts':1}

phys_lb = [-1.0, -1.8, -1.8 ]
phys_ub = [4.2, 2.5, 3.5]
samp_lb = [-6.55, -1.2]
samp_ub = [-6.55, -1.2]
gridsize = [1,1]



result_strings = []
for i in range(n_datasets):
    fitmodel = monod.cme_toolbox.CMEModel('Bursty','None')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,                dataset_strings[i],fitmodel,use_lengths = False,
                gradient_params = gradient_params)
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names)
    full_result_string = inference_parameters.fit_all_grid_points(50,search_data)
    result_strings.append(full_result_string)




gridsize = [6,7]

result_strings = []
for i in range(n_datasets):
    fitmodel = monod.cme_toolbox.CMEModel('Bursty','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = gradient_params)
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names,)
    full_result_string = inference_parameters.fit_all_grid_points(50,search_data)
    result_strings.append(full_result_string)


phys_lb = [-1.8, -1.8,]  # [-7.5,-5] CHECK LENGTH BIAS !!!!! 
phys_ub = [2.5, 2.5] # another check !!!!
samp_lb = [-6.55, -1.2] 
samp_ub = [-6.55, -1.2] 
gridsize = [6,7]

result_strings = []
for i in range(n_datasets):
    fitmodel = monod.cme_toolbox.CMEModel('Constitutive','None')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,                dataset_strings[i],fitmodel,use_lengths = False,
                gradient_params = gradient_params)
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names,)
    full_result_string = inference_parameters.fit_all_grid_points(50,search_data)
    result_strings.append(full_result_string)
    
    
result_strings = []
for i in range(n_datasets):
    fitmodel = monod.cme_toolbox.CMEModel('Constitutive','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = gradient_params)
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names, )
    full_result_string = inference_parameters.fit_all_grid_points(50,search_data)
    result_strings.append(full_result_string)
