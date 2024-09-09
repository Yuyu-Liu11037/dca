import anndata as ad
import pandas as pd
import scipy.sparse
import sys
from tqdm import tqdm

SITE1_CELL = 16311
SITE2_CELL = 25171
SITE3_CELL = 32029
SITE4_CELL = 16750
GEX = 13953

adata = ad.read_h5ad('/workspace/dca/data/citeseq_processed.h5ad')
adata.var_names_make_unique()

X = adata.X.toarray()
X[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX] = 0
X = X.T

# Separate gene and cell information
gene_names = pd.DataFrame(adata.var_names, columns=['Gene'])
cell_names = pd.DataFrame(adata.obs_names, columns=['Cell'])

# Save gene and cell information to separate CSV files
gene_names.to_csv('/workspace/dca/data/citeseq_genes.csv', index=False)
cell_names.to_csv('/workspace/dca/data/citeseq_cells.csv', index=False)

# Save data into chunks
chunk_size = 1000

with tqdm(total=X.shape[0]) as pbar:
    for start in range(0, X.shape[0], chunk_size):
        data_chunk = pd.DataFrame(X[start:start+chunk_size, :], index=gene_names.iloc[start:start+chunk_size, 0], columns=adata.obs_names)
        data_chunk.to_csv('/workspace/dca/data/citeseq_data.csv', mode='a', header=(start == 0))
        pbar.update(chunk_size)

