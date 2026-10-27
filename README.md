# Gut Microbiome Topological Data Analysis Pipeline

This repository contains the analysis code associated with the preprint *Evolution of novel host-microbiota coupling underlies adaptation to nutrient-scarce environments* (bioRxiv, 2026; https://doi.org/10.64898/2026.10.10.758241).  
Code archive (Zenodo): https://doi.org/10.5281/zenodo.23680629
The pipeline applies Mapper-based topological data analysis (TDA) and sparse modeling to identify associations between gut microbiome (or pathway) composition and host morphometric traits.

---

## Scripts

| File | Language | Description |
|------|----------|-------------|
| `DiffDiet_1_mapper_graph.py` | Python | Generates a Mapper graph from ASV abundance data and a morphometric trait lens |
| `DiffDiet_2_mapper_shortestPath.py` | Python | Finds the shortest path between two user-specified nodes in a Mapper graph |
| `DiffDiet_3_elasticnet_mapper.ipynb` | Python | ElasticNetCV sparse modeling, with optional Mapper node-based sample filtering |
| `DiffDiet_3_lassoglm_mapper.mlx` | MATLAB | LassoGLM sparse modeling, with optional Mapper node-based sample filtering |

---

## Provided Data

The code is demonstrated using the **DiffDiet** dataset (`ASV_DiffDiet.csv`), in which fish were raised on different diets (BrineShrimp, Spirulina, Guano).

The same pipeline was also applied to **Antibiotics** and **Fasting** experiment datasets in the paper. To apply these scripts to those datasets, update the parameters as described in the **Adapting to Other Datasets** section below.

---

## Note on the Two Sparse Modeling Scripts

Both DiffDiet_3_elasticnet_mapper.ipynb (Python/ElasticNetCV) and DiffDiet_3_lassoglm_mapper.mlx (MATLAB/lassoglm) use sparse regression to find features linked to a morphometric trait. However, they rely on different algorithms and implementations.

Because of differences in optimization methods and cross-validation procedures, the results (such as selected features and coefficients) may not be exactly the same. Both analyses were carried out and reported in the paper.

---

## Requirements

### Python
```
numpy
pandas
scikit-learn
kmapper
dyneusr
seaborn
matplotlib
scipy
statsmodels
tkinter
dill
```

Install with:
```bash
pip install numpy pandas scikit-learn scikit-tda kmapper dyneusr seaborn matplotlib scipy statsmodels dill
```

### MATLAB
- MATLAB R2021a or later recommended
- Parallel Computing Toolbox (required for `parfor` in `DiffDiet_3_lassoglm_mapper.mlx`; replace `parfor` with `for` if unavailable)
- Python 3.11 environment accessible from MATLAB (required for Mapper node filtering in section 0-2)

---

## Usage

### Step 1 — Generate Mapper graph (`DiffDiet_1_mapper_graph.py`)

Edit the **USER-DEFINED PARAMETERS** section at the top of the script:

```python
INPUT_CSV    = "ASV_DiffDiet.csv"   # Input CSV filename (must be in the same folder as the script)
pop          = "Mo"                  # Population: "SF", "Pa", "Mo", "F1sp", "F1sm", or "all"
morph_trait  = "std_MidgutWidth"     # Morphometric trait used as the Mapper lens
n_cubes      = 4                     # Number of cover intervals
```

Run the script. Outputs are saved in the same folder:
- `dyneusr_*.html` — interactive Mapper graph (DyNeuSR)
- `km_*.html` — Mapper graph (KeplerMapper)
- `nodes_*.pickle` — node membership data (used in subsequent steps)

---

### Step 2 — Find shortest path (`DiffDiet_2_mapper_shortestPath.py`)

Edit the **USER-DEFINED PARAMETERS** section:

```python
INPUT_CSV = "ASV_DiffDiet.csv"   # Same CSV used in Step 1
node1     = "cube0_cluster1"     # Start node (check the Mapper HTML output for node names)
node2     = "cube3_cluster1"     # End node
```

Run the script. A file dialog will open to select the `nodes_*.pickle` file from Step 1.  
The shortest path is saved as a `.txt` file in the same folder as the pickle file.

---

### Step 3a — ElasticNetCV sparse modeling (`DiffDiet_3_elasticnet_mapper.ipynb`)

Open the notebook in Jupyter and edit **Cell 1 (User-Defined Parameters)**:

```python
INPUT_CSV       = "ASV_DiffDiet.csv"
pop2            = "SF"
morph_trait     = "std_HindgutLength"
gut             = "intestine"
LOO             = False        # True = Leave-One-Out CV; False = 10-fold CV × 10 repeats
consistencyRate = 0.8
```

Run cells sequentially. After inspecting the visualization outputs in **Cell 5**, set `l1_selected` in **Cell 6** and continue.  
Optionally, a Mapper nodes pickle file and a cluster info `.txt` file can be provided in **Cell 3** to restrict analysis to a specific subset of Mapper nodes.

---

### Step 3b — LassoGLM sparse modeling (`DiffDiet_3_lassoglm_mapper.mlx`)

Open the Live Script in MATLAB and set MATLAB's **Current Folder** to the folder containing the Live Script and `ASV_DiffDiet.csv` (the script reads the CSV from the current folder, `pwd`). Then edit the **USER-DEFINED PARAMETERS** block in section **0-1**:

```matlab
file_name     = 'ASV_DiffDiet.csv';   % Input CSV filename (must be in MATLAB's Current Folder)
gut           = 'intestine';
pop2          = 'SF';                  % Population: 'SF'  'Pa'  'Mo'  'F1sp'  'F1sm'
s             = 'StandardLength';      % Morphometric trait (response variable Y)
appearanceRate = 0.8;                  % Consistency rate threshold
```

Run sections sequentially using **Ctrl+Enter (Run Section)**. Do not use "Run All" — intermediate plots must be inspected before proceeding.

Key decision point: after running sections **1.1-1** and **1.1-2**, inspect the plots, then set `Selected_L1rate` in section **1.1-3** before continuing.

Optionally, a Mapper nodes pickle file and a cluster info `.txt` file can be provided in section **0-2** to restrict the analysis to a specific subset of Mapper nodes. To use all samples, run section 0-2 and click **Cancel** in the file dialog (or skip section 0-2).


---

## Adapting to Other Datasets

The code can be applied to other datasets (e.g., Antibiotics, Fasting) by modifying the following. The same changes apply to both the Python and MATLAB scripts.

**1. Input file**  
Replace the input filename with your dataset.

**2. Last metadata column**  
The following line marks the boundary between metadata columns and feature columns (ASV or pathway). Update the column name to match the **last metadata column** in your dataset (i.e., the column immediately before the first feature column).

Python:
```python
morp_ind = df_all.columns.get_loc('std_MidgutWidth')
```
MATLAB:
```matlab
morp_ind = find(strcmp(T.Properties.VariableNames, 'std_MidgutWidth'));
```

**3. Population and morphometric trait**  
Update `pop` / `pop2` and `morph_trait` / `s` to values present in your dataset.

**4. Diet labels (scatter plot sections)**  
The scatter plots use diet labels specific to the DiffDiet experiment. Update these to match the condition labels in your dataset.

Python (`DiffDiet_3_elasticnet_mapper.ipynb`, Cell 7):
```python
BS_ind = diet == 'BrineShrimp'
SP_ind = diet == 'Spliru'
GN_ind = diet == 'Guano'
```
MATLAB (`DiffDiet_3_lassoglm_mapper.mlx`, section 1.1-5):
```matlab
ind_BS = ismember(T.diet, 'BrineShrimp');
ind_SP = ismember(T.diet, 'Spliru');
ind_GU = ismember(T.diet, 'Guano');
```

---

## Applying to Pathway Data

The same pipeline can be applied to pathway abundance data (e.g., pathways inferred from ASVs using PICRUSt) with minimal changes. The input CSV has the same format as the ASV table, with pathway columns replacing ASV columns. The only changes required are:

**1. Input file**: replace with your pathway CSV filename.

**2. Last metadata column**: update `'std_MidgutWidth'` to the last metadata column in the pathway table (same position as in the ASV table if the metadata columns are identical).

**3. Column name references**: the scripts use `"ASV_00001"` to verify that feature columns are present. Replace this with the name of the first pathway column in your dataset.

Python (`DiffDiet_1_mapper_graph.py` and `DiffDiet_3_elasticnet_mapper.ipynb`):
```python
df = df.dropna(subset=["ASV_00001"])  # replace with first pathway column name
```

No other changes are needed. All preprocessing (log2 transformation, removal of features present in only one sample) applies equally to pathway data.

---

## Output Files

| File | Description |
|------|-------------|
| `dyneusr_*.html` | Interactive Mapper graph |
| `km_*.html` | KeplerMapper Mapper graph |
| `nodes_*.pickle` | Node membership indices |
| `*_minPath_*.txt` | Shortest path node sequence |
| `ENCV_*_ASVs.csv` | Feature abundance data for selected samples (Python) |
| `ENCV_*_Coeffs.csv` | ElasticNet coefficients and statistics per feature (Python) |
| `ENCV_*_results.pkl` | Full session dump (reloadable with `dill`) (Python) |
| `ML_*_ASVs.csv` | Feature abundance data for selected samples (MATLAB) |
| `ML_*_coeffs.csv` | LassoGLM coefficients and statistics per feature (MATLAB) |
| `ML_*.mat` | Full workspace dump (MATLAB) |

---

## Citation

If you use this code, please cite:  
Evolution of novel host-microbiota coupling underlies adaptation to nutrient-scarce environments. *bioRxiv* (2026). https://doi.org/10.64898/2026.10.10.758241
