

import pickle, os, numbers

import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
from preprocess import read_dataset, normalize

def normalize_tr(x,y):
    adata = sc.AnnData(x)
    adata = read_dataset(adata, transpose=False, test_split=False,copy=True)
    adata= normalize(adata,filter_min_counts_genes=True,size_factors=True, logtrans_input=True)
    sf= adata.obs['size_factors']
    obsnamestr=np.array(adata.obs_names).astype(int)
    y_train=[y[i] for  i in  obsnamestr]
    return sf,adata.raw.X,adata.X,adata.var_names,y_train

def normalize_te(x,var_names,y):
    var_names=np.array(var_names).astype(int)
    x=x[:,var_names]
    adata = sc.AnnData(x)
    adata = read_dataset(adata, transpose=False, test_split=False,copy=True)
    adata= normalize(adata,filter_min_counts_genes=False,size_factors=True, logtrans_input=True)
    sf= adata.obs['size_factors']
    obsnameste=np.array(adata.obs_names).astype(int)
    y_test=[y[i] for  i in  obsnameste]
    
    return sf,adata.raw.X,adata.X,y_test