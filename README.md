# 🧬 Protein Complexity Study Toolkit

This project is part of my Master's thesis in Computer Science. It analyzes the algorithmic complexity of amino acid sequences that form known proteins, with the goal of identifying a simple complexity function that can distinguish between functional protein sequences and non-functional random sequences.


## 📦 Requirements

### Python (>= 3.11 recommended)

Install Python dependencies:

```bash
pip install -r requirements.txt
```


### R Dependencies

Required only for running the Online Algorithmic Complexity Calculator metrics implemented in R.

Ensure you have R (>= 4.5.0) installed.

To install the required R package:

```R
install.packages("acss")
```


## ▶️ Usage

This project can work with any FASTA dataset of amioacid sequences. It can also be adapted for nucleic acids with minimal changes¹.

We used Reviewed (Swiss-Prot) from [UniProt](https://www.uniprot.org/help/downloads).

- Example use cases can be found in `workspace.py`.
- Full function documentation is available in the appendix of `Tesis_Mariano_Oca.pdf`.


### ⚙️ Running Recomendations

This project can be run directly from the console. However, many complexity metrics are computationally intensiveand may occupy the terminal for extended periods. For convenience and background execution, we recommend creating a `execute.sh` file like the following:

```bash
#!/bin/bash
python3 workspace.py
```

Then execute it with:

```bash
nohup ./execute.sh &
```


> \[1] Specifically: `discrepancy()` in `complexity_metrics.py` (line 84) and `random_seq()` in `fasta_utils.py` (line 87).


## 📁 Repository Structure

- **`data/`** – Contains all working files in `.fasta` format along with their corresponding sequence length references in `.txt` format. For example, for the dataset `usp_f.fasta`, the file `sizes_usp_f.txt` contains the lengths of each sequence in order of appearance.

- **`OACC-master/`** – Includes part of the [Online Algorithmic Complexity Calculator](http://www.complexity-calculator.com/) project, used to compute Kolmogorov, Bennett, Shannon complexities, second-order entropy and compression length.

- **`results/`** – Stores the output results from all experiments.

- **`complexity_metrics.py`** – Definition of complexity metrics and the `ComplexitySelector` class, which allows selecting a specific metric and retrieving related configuration and metadata.

- **`fasta_utils.py`** – Utilities for easier handling and processing of FASTA files.

- **`main.py`** – Main script that runs the core experimental workflow.

- **`misc_ploting.py`** – Auxiliary plotting functions used for generating visualizations of experiment results.

- **`misc_utils.py`** – Helper functions used internally by `fasta_utils.py`.

- **`polting_boxplot.ipynb`** – Boxplots, mainly visualizing dataset sequence length distributions.

- **`polting.ipynb`** – Pre-defined plots for each complexity metric.

- **`workspace.py`** – Example Python script for interactive experimentation with the tool.
