# import packages
import os,sys
sys.path.insert(0, '../monod/src/')

import numpy as np
import loompy as lp
import pandas as pd
import pickle

# timing
import time

# import Monod
import monod
from monod import preprocess, extract_data, cme_toolbox, inference, analysis


num_cores = 20

dataset_meta = ['B08']
# 
subcluster_names = ['L2/3 IT','L5 IT','L6 IT','L5/6 NP', 'L6 CT', 'L6b']
subcluster_names = [x.replace(' ','').replace('/','') for x in subcluster_names]
dataset_names = ['allen_'+dataset_meta[0]+'_'+y  for y in subcluster_names]  
raw_data_locations = [dataset_meta[0] for y in subcluster_names] 
transcriptome_filepath = './transcriptome/gg_200524_mouse_genome_polyA_cum_1'

spliced_layer = 'spliced'
unspliced_layer = 'unspliced'
gene_attr = 'gene_name'
cell_attr = 'barcode'

attribute_names=[(unspliced_layer,spliced_layer),gene_attr,cell_attr]

loom_filepaths = ['./data/allen/allen_'+x+'_raw.loom' for x in raw_data_locations] 
n_datasets = len(loom_filepaths)

# filtering
allen_membership = pd.read_csv('./data/allen/cluster.membership.csv',skiprows = 1, names=['barcode','cluster_id'])
allen_annot = pd.read_csv('./data/allen/cluster.annotation.csv')
allen_membership['cell_barcode'] = allen_membership['barcode'].str[:16]
allen_membership['sample'] = allen_membership['barcode'].str[-3:]
allen_membership['cluster_id'] = allen_membership['cluster_id'].astype("category")
allen_annot.set_index('cluster_id',inplace=True)
allen_annot_bc = allen_annot.loc[allen_membership['cluster_id']][['cluster_label','subclass_label','class_label']].set_index(allen_membership.index)
meta = pd.concat((allen_membership,allen_annot_bc),axis=1)
omitted_subtypes = ('L6 IT Car3','L5 ET')

cf = []
thr_lb = [1e4]*4


print(len(loom_filepaths),'len loom file paths')
for k in range(len(dataset_meta)):
    print(k)
    filename = loom_filepaths[len(subcluster_names)*k]
    dataset_name = raw_data_locations[len(subcluster_names)*k]
    
    with lp.connect(filename,mode='r') as ds:
        S = ds.layers[spliced_layer][:]
        U = ds.layers[unspliced_layer][:]
        gene_names = ds.ra[gene_attr]
        bcs = ds.ca[cell_attr]
        n_cells = S.shape[1]
        monod.preprocess.knee_plot(S+U,ax1=None,viz=False,thr=thr_lb[k])
        cf_ = ((S+U).sum(0)>thr_lb[k])
        
        n_annot_bcs = (meta['sample']==dataset_name).sum()
        annot_bcs_in_loom = meta[(meta['sample']==dataset_name)]['cell_barcode'].isin(bcs).sum()
        annot_bcs_in_filt_loom = meta[(meta['sample']==dataset_name)]['cell_barcode'].isin(bcs[cf_]).sum()
        print(f'Dataset {dataset_name}. \n\t{len(bcs)} barcodes in loom, {cf_.sum()} pass filter. {n_annot_bcs} in annotations; of these, {annot_bcs_in_loom} in loom and {annot_bcs_in_filt_loom} in filtered loom.')
        if k==0:
            for subcluster in subcluster_names:
                annot_bcs = meta[(meta['sample']==dataset_name) \
                                           & (meta['subclass_label'].str.replace(' ','').str.replace('/','')==subcluster) \
                                           & ~(meta['subclass_label'].isin(omitted_subtypes))]['cell_barcode']
                cf.append(np.isin(bcs,annot_bcs) & cf_)
                print(f'\t{subcluster}: {len(annot_bcs)} cells in annotations. {np.isin(bcs,annot_bcs).sum()} in loom. {cf[-1].sum()} pass filter.')
                
            for subcluster in omitted_subtypes:
                annot_bcs = meta[(meta['sample']==dataset_name) \
                                           & (meta['subclass_label']==subcluster) ]['cell_barcode']
                CF_ = np.isin(bcs,annot_bcs) & cf_
                print(f'\tOmitted -- {subcluster}: {len(annot_bcs)} cells in annotations. {np.isin(bcs,annot_bcs).sum()} in loom. {CF_.sum()} pass filter.')
        


# preprocessing
dir_string,dataset_strings = monod.preprocess.construct_batch(loom_filepaths, \
                                             transcriptome_filepath, \
                                             dataset_names, \
                                             attribute_names=attribute_names,\
                                             batch_location='./fits',meta='allen_B08',batch_id=1,\
                                             n_genes=3000,exp_filter_threshold=1,cf=cf,creator='mc',datestring='241206')

# set up dictionary for results
timing_dictionary = {d : {} for d in dataset_names}

for d in dataset_names:
    timing_dictionary[d]['models'] = ['Bursty',
                                      'Constitutive',
                                      'Extrinsic',
                                      'Delay',
                                      'DelayedSplicing',
                                      ]
    timing_dictionary[d]['runtime'] = []
    timing_dictionary[d]['n_cells'] = []
    timing_dictionary[d]['n_genes'] = []
    
print(timing_dictionary.keys())

# run inference for different models

# technical noise sampling for ALL models
max_iterations=15
samp_lb = [-9, -4] #-7.5, -2
samp_ub = [-4, 1.5] #-5.5, 0
gridsize = [10,11] #run across grid of technical/sampling parameters

# Bursty
print('working on Bursty')
phys_lb = [-2.0, -1.8, -1.8 ] # b, beta, gamma 
phys_ub = [4.2, 2.5, 2.5] 

result_strings = []
for i,d in enumerate(dataset_names):
    t1_ = time.time()
    fitmodel = monod.cme_toolbox.CMEModel('Bursty','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = {'max_iterations':max_iterations,'init_pattern':'moments','num_restarts':2})
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=cf[i])
    full_result_string = inference_parameters.fit_all_grid_points(num_cores,search_data)

    result_strings.append(full_result_string)
    t2_ = time.time()
    
    timing_dictionary[d]['runtime'].append(t2_-t1_)
    timing_dictionary[d]['n_cells'].append(search_data.n_cells)
    timing_dictionary[d]['n_genes'].append(search_data.n_genes)
    
    
# Constitutive
print('working on Constitutive')
phys_lb = [-1.8, -1.8 ] # beta, gamma
phys_ub = [2.5, 2.5] 

result_strings = []
for i,d in enumerate(dataset_names):
    t1_ = time.time()
    fitmodel = monod.cme_toolbox.CMEModel('Constitutive','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = {'max_iterations':max_iterations,'init_pattern':'moments','num_restarts':2})
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=cf[i])
    full_result_string = inference_parameters.fit_all_grid_points(num_cores,search_data)

    result_strings.append(full_result_string)
    t2_ = time.time()
    
    timing_dictionary[d]['runtime'].append(t2_-t1_)
    timing_dictionary[d]['n_cells'].append(search_data.n_cells)
    timing_dictionary[d]['n_genes'].append(search_data.n_genes)
    
# Extrinsic
print('working on Extrinsic')
phys_lb = [-1.8, -1.8, -1.8 ] # alpha, beta, gamma
phys_ub = [2.5, 2.5, 2.5] 

result_strings = []
for i,d in enumerate(dataset_names):
    t1_ = time.time()
    fitmodel = monod.cme_toolbox.CMEModel('Extrinsic','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = {'max_iterations':max_iterations,'init_pattern':'moments','num_restarts':2})
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=cf[i])
    full_result_string = inference_parameters.fit_all_grid_points(num_cores,search_data)

    result_strings.append(full_result_string)
    t2_ = time.time()
    
    timing_dictionary[d]['runtime'].append(t2_-t1_)
    timing_dictionary[d]['n_cells'].append(search_data.n_cells)
    timing_dictionary[d]['n_genes'].append(search_data.n_genes)
    
    
# Delay
print('working on Delay')
phys_lb = [-1.0, -1.8, -2.5 ] # b, beta, 1/tau
phys_ub = [4.2, 2.5, 2.5] 

result_strings = []
for i,d in enumerate(dataset_names):
    t1_ = time.time()
    fitmodel = monod.cme_toolbox.CMEModel('Delay','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = {'max_iterations':max_iterations,'init_pattern':'moments','num_restarts':2})
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=cf[i])
    full_result_string = inference_parameters.fit_all_grid_points(num_cores,search_data)

    result_strings.append(full_result_string)
    t2_ = time.time()
    
    timing_dictionary[d]['runtime'].append(t2_-t1_)
    timing_dictionary[d]['n_cells'].append(search_data.n_cells)
    timing_dictionary[d]['n_genes'].append(search_data.n_genes)
    
    
# DelayedSplicing
print('working on DelayedSplicing')
phys_lb = [-1.0, -2.5, -1.8 ] # b, 1/tau, gamma
phys_ub = [4.2, 2.5, 2.5] 

result_strings = []
for i,d in enumerate(dataset_names):
    t1_ = time.time()
    fitmodel = monod.cme_toolbox.CMEModel('DelayedSplicing','Poisson')
    inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
                dataset_strings[i],fitmodel,use_lengths = True,
                gradient_params = {'max_iterations':max_iterations,'init_pattern':'moments','num_restarts':2})
    search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
                dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=cf[i])
    full_result_string = inference_parameters.fit_all_grid_points(num_cores,search_data)

    result_strings.append(full_result_string)
    t2_ = time.time()
    
    timing_dictionary[d]['runtime'].append(t2_-t1_)
    timing_dictionary[d]['n_cells'].append(search_data.n_cells)
    timing_dictionary[d]['n_genes'].append(search_data.n_genes)
    
    
    
    
# save timing dictionary
with open('results/allen_B08_timing.pickle', 'wb') as handle:
    pickle.dump(timing_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

