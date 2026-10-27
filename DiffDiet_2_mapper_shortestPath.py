"""
Shortest Path Finder for Mapper Graph
======================================
This script finds the shortest (minimum-weight) path between two user-specified
nodes in a Mapper graph produced by DiffDiet_1_mapper_graph.py.

Workflow:
    1. The ASV CSV file is loaded from the same folder as this script.
    2. A file dialog opens to select the Mapper nodes pickle file.
    3. If population and morphometric trait cannot be inferred from the pickle filename,
       dialog boxes prompt the user for that information.
    4. The shortest path between node1 and node2 is computed and saved as a .txt file.

Node name format:
    Nodes are named by their cover interval and cluster index, e.g.:
        "cube0_cluster1"  ->  interval 0, cluster 1
        "cube3_cluster1"  ->  interval 3, cluster 1

Output:
    A .txt file listing the node sequence of the shortest path,
    saved in the same directory as the input pickle file.
"""

import numpy as np
import pandas as pd
import pickle
import tkinter as tk
import tkinter.filedialog
import tkinter.simpledialog
import os
import itertools as IT

# =============================================================================
# USER-DEFINED PARAMETERS
# =============================================================================

# The two Mapper nodes to find the shortest path between.
# Node name format: "cube{interval_index}_cluster{cluster_index}"
# Check the Mapper HTML output to identify the node names of interest.
node1 = "cube0_cluster1"
node2 = "cube3_cluster1"

# Filename of the ASV CSV file (must be in the same folder as this script).
INPUT_CSV = "ASV_DiffDiet.csv"

# =============================================================================
# DATA LOADING
# =============================================================================

# Build the CSV path relative to this script's location.
# os.path.join avoids issues with Windows backslash path separators.
src = os.path.dirname(os.path.abspath(__file__))
df_all = pd.read_csv(os.path.join(src, INPUT_CSV), index_col=False, header=0)
df = df_all.query('gut == "intestine"')

# =============================================================================
# LOAD MAPPER NODES PICKLE FILE (file dialog)
# =============================================================================

print("Please select the Mapper nodes pickle file.")
tk.Tk().withdraw()
args = {"filetypes": [("pickle", "*.pickle")], "title": "Please select nodes file of pickle format."}
f_path = tkinter.filedialog.askopenfilename(**args)
if len(f_path) == 0:
    raise SystemExit("No pickle file selected. Stopping.")
with open(f_path, mode='br') as fi:
    node_to_csv0 = pickle.load(fi)
file = os.path.splitext(os.path.basename(f_path))[0]

# =============================================================================
# INFER OR PROMPT FOR MORPHOMETRIC TRAIT
# =============================================================================

# Try to detect the morphometric trait from the pickle filename
morp_ind = df.columns.get_loc('std_MidgutWidth')
morp = df.columns[:morp_ind + 1]
flag_m = 0
for i in morp:
    if (i in file):
        flag_m = 1
        morp2 = i

if flag_m != 1:
    # Could not detect from filename -> ask the user
    tk.Tk().withdraw()
    morp2 = tk.simpledialog.askstring('input morph of interest', 'input morph of interest. if not, cancel this')

if type(morp2) is str and len(morp2) > 0:
    df = df.dropna(subset=[morp2])

# =============================================================================
# INFER OR PROMPT FOR POPULATION
# =============================================================================

# Try to detect the population from the pickle filename
pop = df['pop2'].unique()
flag_p = 0
for i in pop:
    if (i in file):
        flag_p = 1
        pop2 = i

if flag_p != 1:
    # Could not detect from filename -> ask the user
    tk.Tk().withdraw()
    pop2 = tk.simpledialog.askstring('input pop of interest', 'input pop of interest. if not, cancel this')

if type(pop2) is str and len(pop2) > 0:
    df = df.query('pop2 == "' + pop2 + '"')
else:
    df = df.copy()

# =============================================================================
# GRAPH CONSTRUCTION FROM MAPPER NODES
# =============================================================================

# Set the morphometric trait values as node weights (mean value per node)
y = df[morp2]

# Reset index to align mapper node member indices with dataframe rows
y2 = y.reset_index()

# Build a DataFrame of nodes and their member sample indices
node_to_csv = list(node_to_csv0.items())
node_df = pd.DataFrame(node_to_csv, columns=['Node_name', 'Members'])
mem_index = tuple(list(node_df['Members']))

# Assign each node a weight = mean morphometric value of its members
node_df['weights'] = [np.mean(y2.iloc[i, 1]) for i in mem_index]

# Find edges: two nodes are connected if they share at least one member sample.
# Edge weight = number of shared members (overlap size).
S = node_df.iloc[:, 1]
member_to_set = [set(ele) for ele in S]
idx = range(len(member_to_set))
pairs = IT.combinations(idx, 2)
nt = lambda a, b: member_to_set[a].intersection(member_to_set[b])
res = dict([(t, nt(*t)) for t in pairs])
res_no_null = {k: v for k, v in res.items() if v}
res_len = {k: len(v) for k, v in res.items() if v}
res_len = {(str(k[0]), str(k[1])): v for k, v in res_len.items()}

# Build adjacency dict indexed by integer-string node indices
graph_mp = {}
for k, v in res_len.items():
    graph_mp.setdefault(k[0], {})[k[1]] = v

# Replace integer-string indices with actual node names (e.g. "cube0_cluster1")
kd = dict(zip(node_df.index.map(str), node_df.Node_name))

def replace_keys(old_dict, key_dict):
    """Recursively replace keys in a nested dict using key_dict mapping."""
    new_dict = {}
    for key in old_dict.keys():
        new_key = key_dict.get(key, key)
        if isinstance(old_dict[key], dict):
            new_dict[new_key] = replace_keys(old_dict[key], key_dict)
        else:
            new_dict[new_key] = old_dict[key]
    return new_dict

nd = replace_keys(graph_mp, kd)

# =============================================================================
# SHORTEST PATH ALGORITHM
# =============================================================================

def find_all_paths(graph, start, end, path=[]):
    """Recursively find all simple paths from start to end in graph."""
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


def min_path(graph, start, end):
    """
    Find the path from start to end whose total edge weight is minimal.
    Edge weight = number of shared members between adjacent nodes.
    Prints all candidate paths and returns the node sequence of the best path.
    """
    paths = find_all_paths(graph, start, end)
    mt = 10**99
    mpath = []
    print('\tAll paths:', paths)
    if len(paths) == 0:
        raise SystemExit(f"No path found between {start} and {end}. "
                         "Check the node names in the Mapper HTML output.")
    for path in paths:
        t = sum(graph[i][j] for i, j in zip(path, path[1::]))
        print('\t\tevaluating:', path, t)
        if t < mt:
            mt = t
            mpath = path

    e1 = ' '.join('{}->{}:{}'.format(i, j, graph[i][j]) for i, j in zip(mpath, mpath[1::]))
    e2 = str(sum(graph[i][j] for i, j in zip(mpath, mpath[1::])))
    print('Best path: ' + e1 + '   Total: ' + e2 + '\n')
    return mpath

# =============================================================================
# RUN AND SAVE RESULTS
# =============================================================================

if __name__ == "__main__":
    graph = nd
    best_path = min_path(graph, node1, node2)

    # Output filename: <pickle name>_minPath_<start node>_<end node>.txt
    # (e.g. "cube0_cluster1" is shortened to "cu0_cl1").
    # If the population or trait was typed in a dialog (not found in the
    # pickle filename), it is added to the output filename, because the
    # sparse modeling scripts check the population/trait in this filename.
    save_name = file
    if flag_p != 1 and type(pop2) is str and len(pop2) > 0:
        save_name += '_' + pop2 + '_'
    if flag_m != 1 and type(morp2) is str and len(morp2) > 0:
        save_name += '_' + morp2 + '_'
    shorten = lambda n: n.replace("be", "").replace("uster", "")
    save_name += '_minPath_' + shorten(best_path[0]) + "_" + shorten(best_path[-1])
    while '__' in save_name:
        save_name = save_name.replace('__', '_')

    save_path = os.path.join(os.path.dirname(f_path), save_name + '.txt')
    pd.Series(best_path).to_csv(save_path, index=False, header=False)
    print(f"Shortest path saved to: {save_path}")
