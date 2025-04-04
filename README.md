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
Two files have not been uploaded to the github due to their size, but are both publically available. The first file is training data from the BindingDB database (https://www.bindingdb.org/rwd/bind/index.jsp). The second missing file is metabolomic data from Campos, A. I. & Zampieri, M. Metabolomics-Driven Exploration of the Chemical Drug Space to Predict Combination Antimicrobial Therapies. Mol. Cell 74, 1291-1303.e6 (2019). Supplementary Table S1 (47MB) can be downloaded and used as input in M2D2. The M2D2 code uses the name "campos.xlsx" for clarity. No preprocessing of the dataset is needed. 
  
