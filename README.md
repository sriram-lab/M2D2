# M2D2

Ultra-high-throughput screening of antimicrobial combination therapies using a two-stage transparent machine learning model M2D2

Here, we present M2D2, a two-stage machine learning (ML) pipeline that identifies promising antimicrobial drug combinations, which are crucial for combating drug resistance. M2D2 addresses key challenges in drug combination discovery by predicting drug synergies using computationally generated drug-protein interaction data, thereby circumventing the need for expensive omics data. 

STAGE 1 ML: generate drug - protein interactions to use as features in stage 2 ML that predicts drug - drug interactions. The two encoding files mine encodings from DeepPurpose library (Huang, K. et al. DeepPurpose: a deep learning library for drug–target interaction prediction. Bioinformatics 1–6 (2020)). To install and access their library of encodings please see their GitHub repository https://github.com/kexinhuang12345/DeepPurpose. Installation instructions found in README.md.
 
Generate target encodings (generate_target_encodings.ipynb) 
	input: .xlsx file containing protein targets and their amino acid sequences
	output: .pkl file containing protein PseudoAAC encodings
	
Generate compound encodings (generate_drug_encodings.ipynb
	input: .xlsx file containing drugs and their SMILES structures
	output: .pkl file containing drug MACCS encodings
	
Generate drug - protein interaction predictions (drug-protein_prediction.ipynb)
	input: bindingdb_merged_all.pkl (training data), drugs_MACCS.pkl (drug encodings), ecoli_4087_pseudoAAC.pkl (protein encodings)
	output: drug - protein interaction values 

STAGE 2 ML: generate drug interaction predictions using any of the available datasets. Training data and stage 1 ML features are already available to use, but can be replaced by user input. The traintest folder provides code to test model performance using different types of data inputs. The predict folder provides all files needed to predict drug interactions for a set of FDA approved drugs combined with training drugs. 

# Illustrations

Provides more extensive heatmaps for all computational and omics datasets. The datasets are protein - drug interactions calculated by ML and molecular docking, and a set of omics data taken from literature including chemogenomics, metabolomics, and transcriptomics. The heatmaps are separated into drugs and pathways. The drugs folder contains individual heatmaps for the 16 drugs that overlap between the five datasets options provided. The pathways folder contains heatmaps for classes of pathways for all drugs and all datasets. 

These visuals are also available as an interactive website: https://sriramlab.shinyapps.io/shiny1/

# NOTE ON LARGE DATASETS
One file has not been uploaded to the github due to their size, but are both publically available. The missing file is metabolomic data from Campos, A. I. & Zampieri, M. Metabolomics-Driven Exploration of the Chemical Drug Space to Predict Combination Antimicrobial Therapies. Mol. Cell 74, 1291-1303.e6 (2019). Supplementary Table S1 (47MB) can be downloaded and used as input in M2D2. The M2D2 code uses the name "campos.xlsx" for clarity. No preprocessing of the dataset is needed. 
  

# HOW TO RUN

## Requirements

Stage 1 uses Python. The main Stage 2 analyses use MATLAB, with a Python notebook provided as a simplified Stage 2 demo.

Install the basic Python dependencies with:

```bash
pip install numpy pandas scikit-learn scipy openpyxl matplotlib
```

Stage 1 also requires:

* [DeepPurpose](https://github.com/kexinhuang12345/DeepPurpose) for PseudoAAC protein encodings
* `rdkit` for MACCS drug encodings

The full Stage 2 workflows in `M2D2_stage2_predict/` and `M2D2_stage2ML_traintest/` require MATLAB with the Statistics and Machine Learning Toolbox.

## Quick start: Stage 2 demo

The easiest way to run M2D2 is the Python Stage 2 demo:

```text
M2D2_demo_python/M2D2_stage2_solution.ipynb
```

Run the notebook from inside `M2D2_demo_python/`. All required inputs are provided in `input_data/`, so this demo does not require MATLAB, DeepPurpose, or rerunning Stage 1.

`M2D2_stage2_demo.ipynb` contains the same workflow with the two main M2D2 functions left blank as an exercise.

A MATLAB version of the demo is available in `M2D2_demo_matlab/`.

## Stage 1: Drug–protein interaction prediction

Stage 1 generates predicted drug–protein interactions for use as features in Stage 2.

Run the notebooks from inside:

```text
M2D2_stage1ML/
```

The workflow consists of:

| Step | Notebook                           | Output                     |
| ---- | ---------------------------------- | -------------------------- |
| 1    | `generate_drug_encodings.ipynb`    | `drugs_MACCS.pkl`          |
| 2    | `generate_protein_encodings.ipynb` | `ecoli_4087_pseudoAAC.pkl` |
| 3    | `drug-protein_prediction.ipynb`    | `out.csv`                  |

The drug and protein encodings from steps 1 and 2 are already included in the repository. To rerun the drug–protein prediction, you can therefore begin with step 3.

The prediction model is trained on BindingDB using MACCS drug fingerprints and PseudoAAC protein encodings, and then predicts binding affinity for each drug–protein pair.

A precomputed Stage 1 feature matrix is already included at:

```text
M2D2_stage2_predict/datasets/ml_4070x58.xlsx
```

Therefore, Stage 2 can be run without rerunning Stage 1.

## Stage 2: Drug–drug interaction prediction

The full Stage 2 workflows are implemented as MATLAB Live Scripts.

Set the MATLAB working directory to the folder containing the workflow you want to run.

### Model evaluation

To evaluate M2D2 using the available feature datasets:

```text
M2D2_stage2ML_traintest/MLAnalysis_weighted.mlx
```

### Drug interaction prediction

To predict interactions involving the FDA-approved drug set:

```text
M2D2_stage2_predict/MLAnalysis_test2000drugList.mlx
```

Stage 2 converts the drug–protein feature matrix into binary interaction profiles, constructs the M2D2 sigma and delta features for each drug pair, and trains a random forest model to predict drug–drug interaction scores.

The default Stage 1 ML features can be replaced with the other feature datasets included in the repository, including chemogenomics, molecular docking, metabolomics, and transcriptomics.

Training drug–drug interaction scores are provided in:

```text
drugInteractions/drugInteractions_weights.xlsx
```

## Boltz-2

`M2D2_Boltz2/` provides an alternative Stage 1 implementation using [Boltz-2](https://github.com/jwohlwend/boltz) for drug–protein predictions. The resulting Boltz-2 affinity scores are used as features for the M2D2 Stage 2 drug–drug interaction model.

A precomputed Boltz-2 feature matrix is included, so the Stage 2 analyses can be run without rerunning the GPU-intensive Boltz-2 predictions.

See `M2D2_Boltz2/README.md` for setup and usage instructions.


## M. tuberculosis

The `M2D2_mtb/` folder contains the *Mycobacterium tuberculosis* implementation of M2D2.

See `M2D2_mtb/README.md` for instructions specific to that workflow.

