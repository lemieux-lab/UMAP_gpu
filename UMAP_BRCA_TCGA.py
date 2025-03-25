from cuml import UMAP as UMAP_GPU
import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cuml import UMAP as UMAP_GPU
from cuml.decomposition import PCA as PCA_GPU
import pdb

def plot_TCGA(X_tr, labels, clrs):
  fig = plt.figure(dpi = 300)
  axis = fig.add_subplot(111)
  for lab in np.unique(labels):
      group = labels == lab
      col = clrs[clrs["labs"] == lab].hexcolor.values[0]
      axis.scatter(X_tr[group, 0], X_tr[group, 1],s=12, c = col,linewidths=0.5, edgecolors = "k")
  plt.tight_layout()

def plot_BRCA(X_tr, labels):
  fig = plt.figure(dpi = 300)
  axis = fig.add_subplot(111)
  for lab in np.unique(labels):
      group = labels == lab
      axis.scatter(X_tr[group, 0], X_tr[group, 1],s=12,linewidths=0.5, label = lab, edgecolors = "k")
  plt.legend()
  plt.tight_layout()


# fname = "BRCA_TCGA_tpm_pam50_Thennavan"
fname = "BRCA_TCGA_pam50_Thennavan"
TCGA_data = h5py.File(f"{fname}.h5", "r")

tpm_data = np.log10(TCGA_data["Data Matrix"][TCGA_data["Gene Type"][:].astype(str) == "protein_coding",:].T + 1)
labels = TCGA_data["PAM50"][:].astype(str)

umap_gpu = UMAP_GPU(n_neighbors=30, min_dist = 0.7, n_components=2)
X_tr_umap = umap_gpu.fit_transform(tpm_data)

plot_BRCA(X_tr_umap, labels)
plt.savefig(f"UMAP_{fname}.svg")
plt.savefig(f"UMAP_{fname}.png")
