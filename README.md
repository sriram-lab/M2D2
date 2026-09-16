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

Stage 1 uses Python. Stage 2 is available in both Python and MATLAB.

Install the basic Python dependencies with:

```bash
pip install numpy pandas scikit-learn scipy openpyxl matplotlib
```

Stage 1 additionally requires:

* [DeepPurpose](https://github.com/kexinhuang12345/DeepPurpose) for PseudoAAC protein encodings
* `rdkit` for MACCS drug encodings

The full MATLAB Stage 2 workflows require MATLAB with the Statistics and Machine Learning Toolbox.

`bindingdb_merged_all.pkl`, used to train the Stage 1 model, is stored with Git LFS. Make sure Git LFS is installed if you plan to rerun Stage 1.

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

The drug and protein encodings from steps 1 and 2 are already included in the repository, so the drug–protein prediction can be rerun starting from step 3.

The prediction model is trained on BindingDB using MACCS drug fingerprints and PseudoAAC protein encodings, and predicts a binding score for each drug–protein pair.

A precomputed Stage 1 feature matrix is also included at:

```text
M2D2_stage2_predict/datasets/ml_4070x58.xlsx
```

Therefore, Stage 2 can be run directly without rerunning Stage 1.

## Stage 2: Drug–drug interaction prediction

Stage 2 uses the drug–protein features from Stage 1 to construct M2D2 sigma and delta features and predict drug–drug interaction scores.

Stage 2 is implemented in both Python and MATLAB.

### Python

The Python implementation can be run through:

```text
M2D2_demo_python/M2D2_stage2_solution.ipynb
```

Run the notebook from inside `M2D2_demo_python/`. All required inputs are provided in `input_data/`, so Stage 1 does not need to be rerun.

`M2D2_stage2_demo.ipynb` contains the same workflow with the two core M2D2 functions left blank as an exercise.

A MATLAB version of the same demo is available in:

```text
M2D2_demo_matlab/
```

The `M2D2_Boltz2/` workflow also provides a Python implementation of Stage 2, using Boltz-2 predictions in place of the original Stage 1 features. See the Boltz-2 section below for details.

### MATLAB

The original Stage 2 workflows are also provided as MATLAB Live Scripts.

Set the MATLAB working directory to the corresponding folder before running the script.

| Goal                                               | Folder                     | Run                               |
| -------------------------------------------------- | -------------------------- | --------------------------------- |
| Evaluate model performance across feature datasets | `M2D2_stage2ML_traintest/` | `MLAnalysis_weighted.mlx`         |
| Predict interactions for the FDA-approved drug set | `M2D2_stage2_predict/`     | `MLAnalysis_test2000drugList.mlx` |

The Stage 2 workflow binarizes the Stage 1 drug–protein feature matrix, constructs sigma and delta profiles for each drug pair, and trains a random forest model to predict drug–drug interaction scores.

In addition to the Stage 1 ML predictions, the MATLAB workflow can use the other feature datasets included in the repository, including chemogenomics, molecular docking, metabolomics, and transcriptomics.

Training drug–drug interaction scores are provided in:

```text
drugInteractions/drugInteractions_weights.xlsx
```

## Boltz-2

`M2D2_Boltz2/` provides an alternative Stage 1 implementation using [Boltz-2](https://github.com/jwohlwend/boltz) for drug–protein predictions.

The resulting Boltz-2 affinity scores are used as input features for the same M2D2 Stage 2 workflow. The folder includes a Python Stage 2 implementation for evaluating these features.

A precomputed Boltz-2 feature matrix is included, so Stage 2 can be run without rerunning the GPU-intensive Boltz-2 predictions.

See `M2D2_Boltz2/README.md` for setup and usage instructions.

## M. tuberculosis

`M2D2_mtb/` applies the M2D2 workflow to *Mycobacterium tuberculosis*, with organism-specific drug and protein inputs.

A precomputed Stage 1 output is included, so its Stage 2 analysis can also be run directly.

See `M2D2_mtb/README.md` for setup and usage instructions.
