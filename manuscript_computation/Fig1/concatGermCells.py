import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import anndata
import loompy
import scipy.io as sio

data_path = '/home/tchari/counts/germcell_splicing/'
meta_path = '/home/tchari/metadata/'
out_path = "/home/tchari/counts/germcell_splicing/loom/"


testS = pd.read_csv(data_path+'GSE136220_raw_counts_matrix_pgcs_no_adrenal.csv',index_col=0)
testU = pd.read_csv(data_path+'GSE136220_raw_counts_matrix_pgcs_no_adrenal_unspliced.csv',index_col=0)
ids = [s1.split(':')[0] for s1 in testS.index]

meta = pd.DataFrame()
meta['identity'] = ids
meta['cell_barcode'] = list(testS.index)
meta.to_csv(meta_path+'germCell_meta.csv',index=None)

print(meta.head())

geneNames = testS.columns.values.tolist() #Check
gene_meta = pd.read_csv(meta_path+'t2g_mouse.txt',delimiter='\t',header=None)
gene_dict = dict(zip(gene_meta[1], gene_meta[2]))
geneNames = [gene_dict[i] if i in list(gene_dict.keys()) else i for i in geneNames]

S = testS.to_numpy()
U = testU.to_numpy()

fname = out_path+'allconds.loom'

print('Making Loom')
#Make loom of U/S
retAdata = anndata.AnnData(
			X=sparse.csr_matrix(S),
			layers={
				'spliced': sparse.csr_matrix(S),
				'unspliced': sparse.csr_matrix(U),
			},
			obs=pd.DataFrame({'barcode': np.array(list(testS.index))},index=np.array(list(testS.index))),
			var=pd.DataFrame({'gene_name': np.array(geneNames)},index=np.array(geneNames))
		)

retAdata.write_loom(fname)

