# M2D2 for *Mycobacterium tuberculosis*

This folder contains the *Mycobacterium tuberculosis* implementation of the M2D2 pipeline.

The overall workflow is the same as the main M2D2 pipeline:

1. **Stage 1:** predict drug–protein interactions across the *M. tuberculosis* proteome.
2. **Stage 2:** use the predicted drug–protein interactions to predict drug–drug interactions.

See the repository root README for a general description of M2D2.

```text
stage1/   Drug–protein prediction
stage2/   Drug–drug interaction prediction
```

Run scripts and notebooks from their respective directories, since file paths are relative to the working directory.

## Requirements

In addition to the packages listed in the root README:

```bash
pip install xgboost seaborn requests
```

## Stage 1: Drug–protein prediction

The Stage 1 workflow predicts binding affinity between the 52 drugs in `drugs_mtb.xlsx` and 3,996 proteins from the *M. tuberculosis* proteome.

The required drug and protein encodings have already been generated and are included in the repository. To run the prediction:

```bash
cd stage1
jupyter notebook drug-protein_prediction.ipynb
```

The notebook uses:

* `drugs_MACCS.pkl` — MACCS encodings for the 52 drugs
* `mtbSequence.pkl` — protein sequence encodings
* `bindingdb_merged_all.pkl` — BindingDB training data

and produces:

```text
out.csv
```

`bindingdb_merged_all.pkl` is approximately 510 MB and is not included in the repository. See the large-dataset instructions in the root README for obtaining it.

A precomputed `out.csv` is also included in `stage2/`, so Stage 1 can be skipped if you only want to run the drug–drug interaction model.

### Regenerating the Stage 1 inputs

The following scripts/notebooks were used to generate the included drug and protein inputs:

| Step | Script/notebook                    | Output             |
| ---- | ---------------------------------- | ------------------ |
| 1    | `getUniqueDrugNames.py`            | `unique_drugs.txt` |
| 2    | `SMILES_scraper.py`                | `drugs_mtb.xlsx`   |
| 3    | `generate_drug_encodings.ipynb`    | `drugs_MACCS.pkl`  |
| 4    | `generate_protein_encodings.ipynb` | `mtbSequence.pkl`  |

These steps only need to be rerun if the drug set or protein set is changed.

### Binding-affinity cutoff

Stage 1 predicts binding affinity in nM, where lower values indicate stronger predicted binding. Stage 2 therefore retains the bottom 20% of predicted affinities for each drug when constructing the binary drug–protein interaction matrix.

This differs from the original E. coli data, where the stored prediction score is oriented in the opposite direction.

## Stage 2: Drug–drug interaction prediction

To run Stage 2:

```bash
cd stage2
jupyter notebook ML_Analysis.ipynb
```

Run the notebook from top to bottom.

The notebook reads:

* `out.csv` — predicted drug–protein affinities from Stage 1
* `drugs_mtb.xlsx` — drug information and SMILES
* an experimental drug–drug interaction dataset

It then constructs the M2D2 sigma/delta features and trains an XGBoost regression model.

The notebook evaluates the model using:

* repeated 70/30 train/test splits
* leave-one-interaction-out validation
* leave-one-drug-out validation
* leave-one-drug-mechanism-class-out validation

### Interaction datasets

The interaction dataset can be selected near the beginning of `ML_Analysis.ipynb`.

Available files include:

| File                            | Description                                        |
| ------------------------------- | -------------------------------------------------- |
| `all_drugs_duplicates.xlsx`     | All measurements, including replicate measurements |
| `all_drugs_avg.xlsx`            | Replicate measurements averaged                    |
| `all_but_rhoads_duplicate.xlsx` | Rhoads dataset excluded; replicates retained       |
| `all_but_rhoads_avg.xlsx`       | Rhoads dataset excluded; replicates averaged       |

The current notebook configuration uses:

```python
drug_interactions = pd.read_excel("all_drugs_duplicates.xlsx")
```

The interaction data were compiled from Ma et al. (2019), Yilancioglu et al. (2019), and the Rhoads dataset.

## Main model settings

The main parameters that may be useful to modify are:

* **Protein-binding percentile cutoff:** 20%
* **Number of repeated holdout splits:** 100
* **XGBoost hyperparameter grid**
* **Drug mechanism groups**

These settings are defined near the beginning of `ML_Analysis.ipynb`.

## Plotting

`data_graphs_before_model.py` and `data_graphs_after_model.py` contain optional plotting code used with `ML_Analysis.ipynb`.

They depend on variables created by the notebook and are not intended to be run as standalone scripts.
