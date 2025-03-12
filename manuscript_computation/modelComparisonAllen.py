# import packages
import sys
sys.path.insert(0, '../monod/src/')

import glob

import numpy as np
import pickle

# import Monod
import monod
from monod import preprocess, extract_data, cme_toolbox, inference, analysis


def get_logL_AIC(search_results,search_data, EPS=1e-20,offs=0):
    ''' Returns the AIC of all sr_arr genes using logL at the sampling optimum.
    
    parameters
    -----------------------------------------------
    search_results : Monod search result class
    search_data: associated Monod search data class
    
    returns
    -----------------------------------------------
    AIC: vector of length sr_arr.n_genes with AIC values according to the model of search_results
    
    '''
    
    hist_type = inference.get_hist_type(search_data)
    logL = np.zeros(search_results.n_genes)
    
    
    for gene_index in range(search_results.n_genes):
        logL[gene_index] = search_results.model.eval_model_logL(
                p=search_results.phys_optimum[gene_index],
                limits=search_data.M[:, gene_index] + offs,
                samp=search_results.regressor_optimum[gene_index],
                data=search_data.hist[gene_index],
                hist_type=hist_type,
                EPS=EPS,
                n_cells = search_results.n_cells
            )
    
    AIC = 2*search_results.sp.n_phys_pars - 2*logL

    
    return(logL,AIC)


result_paths = [g for g in glob.glob('./fits/mc_240521_026_allen_B08_1/*') if '/allen_B08' in g]
models = ['Bursty','Constitutive','Extrinsic','Delay','DelayedSplicing']



num_cores = 12
# for path in result_paths:
#     save_dict = {}
#     sample = path.split('/')[-1]
#     print(sample)
    
    
#     sd_ = monod.analysis.load_search_data(f'{path}/raw.sd')
#     gene_names_ = sd_.gene_names
    
#     for model in models:
#         print(model)
#         save_dict[model] = {}
#         if model == 'DelayedSplicing':
#             sr_ = monod.analysis.load_search_results(
#                 f'./fits/mc_240521_029_allen_B08_1/{sample}/{model}_Poisson_10x11/grid_scan_results.res')
#         else:
#             sr_ = monod.analysis.load_search_results(f'{path}/{model}_Poisson_10x11/grid_scan_results.res')
#         sr_.find_sampling_optimum()
#         _ = sr_.chisquare_testing(sd_)
#         sr_.chisq_best_param_correction(sd_,viz=False) 
#         sr_.compute_sigma(sd_,num_cores=num_cores)
#         rejected_genes_ = sr_.rejected_genes
#         logL_,AIC_ = get_logL_AIC(sr_,sd_)
        
#         save_dict[model]['logL'] = logL_
#         save_dict[model]['AIC'] = AIC_
#         save_dict[model]['rejected_genes'] = rejected_genes_
#         save_dict[model]['gene_names'] = gene_names_
    
#     with open(f'results/{sample}.pickle', 'wb') as handle:
#         pickle.dump(save_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
models = ['DelayedSplicing']
for path in result_paths:
    sample = path.split('/')[-1]
    with open(f"./results/{sample}.pickle", "rb") as file:
        save_dict = pickle.load(file)

    print(sample)
    
   
    sd_ = monod.analysis.load_search_data(f'{path}/raw.sd')
    gene_names_ = sd_.gene_names
    
    for model in models:
        print(model)
        save_dict[model] = {}
        if model == 'DelayedSplicing':
            sr_ = monod.analysis.load_search_results(
                f'./fits/mc_240521_029_allen_B08_1/{sample}/{model}_Poisson_10x11/grid_scan_results.res')
        else:
            sr_ = monod.analysis.load_search_results(f'{path}/{model}_Poisson_10x11/grid_scan_results.res')
        sr_.find_sampling_optimum()
        _ = sr_.chisquare_testing(sd_)
        sr_.chisq_best_param_correction(sd_,viz=False) 
        sr_.compute_sigma(sd_,num_cores=num_cores)
        rejected_genes_ = sr_.rejected_genes
        logL_,AIC_ = get_logL_AIC(sr_,sd_)
        
        save_dict[model]['logL'] = logL_
        save_dict[model]['AIC'] = AIC_
        save_dict[model]['rejected_genes'] = rejected_genes_
        save_dict[model]['gene_names'] = gene_names_
    
    with open(f'results/{sample}.pickle', 'wb') as handle:
        pickle.dump(save_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)        
        
