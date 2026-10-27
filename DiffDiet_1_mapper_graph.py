"""
Topological Data Analysis (TDA) using Mapper Algorithm
=======================================================
This script generates a Mapper graph from gut microbiome (ASV) data
combined with host morphometric traits.
 
Output:
    - DyNeuSR interactive HTML graph
    - KeplerMapper HTML graph
    - Pickle file of node information
 
Reference:
    KeplerMapper: https://kepler-mapper.scikit-tda.org/
    DyNeuSR:      https://braindynamicslab.github.io/dyneusr/
 
Requirements:
    numpy, pandas, scikit-learn, kmapper, dyneusr, seaborn, matplotlib
"""
 
import os
import pickle
 
import numpy as np
import pandas as pd
import kmapper as km
import sklearn.cluster
from kmapper import KeplerMapper, Cover
from dyneusr import DyNeuGraph
 
# =============================================================================
# USER-DEFINED PARAMETERS
# =============================================================================
 
# Path to the input CSV file (relative to this script, or absolute path)
INPUT_CSV = "ASV_DiffDiet.csv"
 
# Population to analyze.
# Options: "SF", "Pa", "Mo", "F1sp", "F1sm", "all"
pop = "Mo"
 
# Morphometric trait to use as the Mapper lens.
# Options: "StandardLength", "std_MidgutWidth", "std_HindgutLength", "std_TotalGutLength"
morph_trait = "std_MidgutWidth"
 
# Number of intervals (cubes) for the Mapper cover
n_cubes = 4
 
# Percent overlap between intervals for the Mapper cover
perc_overlap = 0.5
 
# Number of clusters for KMeans (used as the Mapper clusterer)
n_clusters = 2
 
# =============================================================================
# DATA LOADING
# =============================================================================
 
file_path = os.path.abspath(__file__)
src_dir = os.path.dirname(file_path)
data_path = os.path.join(src_dir, INPUT_CSV)
 
df_all = pd.read_csv(data_path, index_col=False, header=0)
 
# =============================================================================
# DATA PREPROCESSING
# =============================================================================
 
# Drop samples missing the target morphometric trait or ASV data
df = df_all.dropna(subset=[morph_trait, "ASV_00001"])
 
# Restrict to intestine samples only
df = df.query('gut == "intestine"')
 
# Filter by population if not analyzing all populations
if pop != "all":
    df = df.query(f'pop2 == "{pop}"')
 
# Ensure no remaining NaN in the lens column after filtering
df = df.dropna(subset=[morph_trait])
 
# Separate morphometric metadata (y) from ASV abundance data (X)
morp_ind = df_all.columns.get_loc("std_MidgutWidth")
X = df.iloc[:, morp_ind + 1:]   # ASV columns
y = df.iloc[:, :morp_ind + 1]   # Metadata + morphometric columns
y = y.drop(columns="Pyloric.caeca")
 
# Log2-transform ASV counts; replace -inf (from log2(0)) with 0
X = np.log2(X)
X = X.replace([np.inf, -np.inf], 0)
 
# =============================================================================
# MAPPER GRAPH CONSTRUCTION
# =============================================================================
 
# Initialize KeplerMapper
mapper = KeplerMapper(verbose=1)
 
# Use the selected morphometric trait directly as the 1-D lens
lens = df[morph_trait]
 
# Build the Mapper graph
#   - Cover: n_cubes intervals with perc_overlap overlap
#   - Clusterer: KMeans with n_clusters
graph = mapper.map(
    lens,
    X,
    clusterer=sklearn.cluster.KMeans(n_clusters),
    cover=km.Cover(n_cubes=n_cubes, perc_overlap=perc_overlap),
)
 
# =============================================================================
# OUTPUT
# =============================================================================
 
# Construct output filename based on parameters
clusterer_name = graph["meta_data"]["clusterer"]
base_name = f"{pop}_{morph_trait}_{clusterer_name}_n_cube-{n_cubes}_"
 
dyneusr_html  = f"dyneusr_{base_name}.html"
km_html       = f"km_{base_name}.html"
nodes_pickle  = f"nodes_{base_name}.pickle"
 
# DyNeuSR interactive visualization
dG = DyNeuGraph(G=graph, y=y)
dG.visualize(dyneusr_html)
 
# KeplerMapper visualization
mapper.visualize(graph, path_html=km_html)
 
# Save node membership information as a pickle file
nodes = dG.cache_.get("G").get("nodes")
with open(nodes_pickle, mode="wb") as f:
    pickle.dump(nodes, f)
 
print(f"\nOutputs saved:\n  {dyneusr_html}\n  {km_html}\n  {nodes_pickle}")

    
    



##save the network stracture##
