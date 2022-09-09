import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


counts_matrix = scipy.io.mmread('./filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes('./filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)



barcodes = np.array(pd.read_csv('./filtered_feature_bc_matrix/barcodes.tsv', header=None, index_col=None))


a = np.array([barcodes[:,0], doublet_scores, predicted_doublets])
data = pd.DataFrame({'barcodes': a[0,:], 'scores':a[1,:], 'prediction':a[2,:]})
data.to_csv('./doublet_prediction.txt', index=False, header=True)







