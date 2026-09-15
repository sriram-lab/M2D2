# Boltz-2 as the Stage 1 feature basis for M2D2

This folder replaces M2D2's original Stage 1 drug–target predictions with [Boltz-2](https://github.com/jwohlwend/boltz), while leaving the Stage 2 M2D2 workflow unchanged.

A precomputed Boltz-2 feature matrix is included, so Stage 2 can be run without rerunning the GPU-intensive Stage 1 predictions.

## Layout

This folder uses three input files from the main M2D2 repository:

```text
M2D2/
├── M2D2_stage1ML/
│   ├── drugs.xlsx                            # drug names, abbreviations, and SMILES
│   └── ecoli_4087_AAseq.xlsx                 # protein IDs and sequences
├── M2D2_stage2ML_traintest/
│   └── drugInteractions/
│       └── drugInteractions_weights.xlsx     # drug–drug interaction scores
└── M2D2_Boltz2/
    ├── README.md
    ├── stage1/
    │   ├── predict.py                        # Boltz-2 drug–protein predictions
    │   └── extract.py                        # prediction outputs -> feature matrices
    ├── stage2/
    │   ├── utils.py                          # features, models, and metrics
    │   └── evaluate.py                       # Stage 2 model evaluation
    └── data/
        ├── affinity_probability_binary.csv   # precomputed Stage 1 feature matrix
        └── leaveOneDrugOut_w75_n25_v2.xlsx   # drug mechanism classes
```

Input paths can be specified through command-line arguments if the files are stored elsewhere.

Running Stage 1 and Stage 2 may additionally create:

```text
M2D2_Boltz2/
├── msa_cache/                    # cached protein MSAs
├── run/                          # Boltz-2 inputs and prediction outputs
│   ├── yamls/
│   └── pred/
└── results/                      # Stage 2 evaluation results
```

These generated files are not included in version control.

## Requirements

```bash
pip install pandas numpy scipy scikit-learn pyyaml xgboost
pip install boltz
```

`boltz` is required only for Stage 1 and currently supports Python 3.10–3.12.

Stage 2 can be run with Random Forest alone using `--models rf`, in which case XGBoost is not required.

Boltz-2 should be installed as a Python dependency rather than copied into this repository. `predict.py` uses the installed `boltz` package to run predictions and generate MSAs.

The results included here were generated with Boltz 2.2.1. Pinning the same version is recommended for exact reproduction.

Boltz model weights and CCD data are downloaded automatically on first use. By default they are stored in `~/.boltz`; the location can also be changed with `$BOLTZ_CACHE` or the `--cache` option.

## Stage 1 — Boltz-2 predictions

Stage 1 requires a GPU. The full drug × protein prediction set is computationally expensive, so parallel execution is recommended when multiple GPUs are available.

A basic run is:

```bash
python stage1/predict.py --fetch-msa \
    --drugs ../M2D2_stage1ML/drugs.xlsx \
    --proteins ../M2D2_stage1ML/ecoli_4087_AAseq.xlsx \
    --interactions ../M2D2_stage2ML_traintest/drugInteractions/drugInteractions_weights.xlsx \
    --out run
```

By default, predictions are generated for all drugs present in the interaction dataset.

To restrict the run to selected drugs:

```bash
--agents AMK RIF
```

This can also be used to add new drugs later without recomputing existing predictions.

Completed drug–protein pairs are skipped automatically, so interrupted runs can be resumed by rerunning the same command.

### MSA generation

With `--fetch-msa`, protein MSAs are generated and cached for reuse across drugs and subsequent runs.

MSA generation and Boltz-2 prediction can overlap: proteins whose MSAs are already available can be processed while the remaining MSAs are still being generated.

The default batch size is 8 proteins per Boltz-2 invocation. This can be changed with `--batch-size`:

* `--batch-size 8` — default
* `--batch-size 1` — process proteins as soon as their MSAs become available
* `--batch-size 0` — wait for all MSAs before starting prediction

`--queue-depth` controls how far MSA generation may run ahead of prediction.

To populate the MSA cache without running Boltz-2 predictions:

```bash
python stage1/predict.py --fetch-msa --yaml-only --out run
```

This is recommended before starting multiple parallel prediction processes, since the same protein MSAs are reused across drugs.

### Splitting Stage 1 across multiple GPUs

The protein panel can be divided into independent shards:

```text
--nshards N
--shard I
```

`--nshards N` divides the protein set into `N` groups, and `--shard I` runs one of those groups. Each shard processes all drugs against its assigned subset of proteins.

For example:

```bash
python stage1/predict.py \
    --shard 0 \
    --nshards 4 \
    --out run
```

Run the corresponding commands with `--shard 0`, `1`, `2`, and `3` on separate GPUs to process four shards in parallel.

Each shard writes to its own output directory, so the processes run independently. If one shard is interrupted, it can be rerun using the same `--shard` and `--nshards` values without restarting the others.

MSAs should generally be generated once before starting multiple shards:

```bash
python stage1/predict.py --fetch-msa --yaml-only --out run
```

Afterward, launch the desired number of prediction shards using the shared MSA cache.

### Collecting Stage 1 results

After all predictions are complete, combine the Boltz-2 outputs with:

```bash
python stage1/extract.py \
    --pred run/pred \
    --results data \
    --drop-incomplete-proteins
```

`extract.py` searches the prediction directory recursively, so the same command works for both sharded and non-sharded runs.

One feature matrix is written for each affinity metric reported by Boltz-2. The Stage 2 analyses in this folder use:

```text
data/affinity_probability_binary.csv
```

`--drop-incomplete-proteins` removes proteins missing a prediction for one or more drugs, ensuring that all drugs share the same protein feature set.

`--merge-into` can be used to add newly predicted drugs to an existing feature matrix.

A precomputed `affinity_probability_binary.csv` is included in the repository, so Stage 2 can be run without rerunning Stage 1.

### Boltz-2 prediction settings

The included Stage 1 results were generated with:

```text
--sampling_steps 100
--sampling_steps_affinity 100
--diffusion_samples_affinity 5
```

These settings are defined in `predict.py`.

Other useful options include:

* `--msa-dir` — location of the MSA cache
* `--cache` — location of Boltz model weights and related data
* `--devices` — number of GPUs used by each process
* `--yaml-only` — generate Boltz-2 inputs without running predictions
* `--batch-size` — number of proteins passed to each Boltz-2 invocation
* `--queue-depth` — limit how far MSA generation runs ahead of prediction

## Stage 2 — Drug–drug interaction models

Stage 2 can be run on a standard CPU machine:

```bash
python stage2/evaluate.py \
    --matrix data/affinity_probability_binary.csv \
    --drugs ../M2D2_stage1ML/drugs.xlsx \
    --interactions ../M2D2_stage2ML_traintest/drugInteractions/drugInteractions_weights.xlsx \
    --classes data/leaveOneDrugOut_w75_n25_v2.xlsx \
    --results results
```

Two experiment groups are run by default. Either can be run individually using `--experiment representations` or `--experiment published`.

### Representation comparison

```bash
--experiment representations
```

This experiment compares five ways of constructing M2D2 features from the Boltz-2 affinity predictions:

* unbinarized features
* global 0.5 cutoff
* per-drug 80th-percentile cutoff
* global-cutoff features with the delta score
* per-drug-cutoff features with the delta score

Each representation is evaluated with Random Forest and XGBoost using:

* 5-fold cross-validation over drug combinations
* leave-one-drug-out evaluation

Results are written to:

```text
representations.csv
representations_ldo_per_drug.csv
```

### Published M2D2 protocol

```bash
--experiment published
```

This evaluates the Boltz-2 features using the M2D2 protocols reported in the original workflow:

* 50 repeated 70/30 train/test splits
* leave-one-drug-out evaluation

The repeated holdout evaluation uses a 100-tree Random Forest configured to match the MATLAB `TreeBagger` setup.

Results are written to:

```text
published_protocols.csv
```

## Additional Stage 2 options

Useful options include:

* `--models rf` — Random Forest only
* `--models xgb` — XGBoost only
* `--percentile` — percentile used for per-drug binarization
* `--iterations` — number of repeated holdout evaluations
* `--seed` — random seed used for cross-validation
