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

This project can work with any FASTA dataset of amioacid sequences. It can also be adapted for nucleic acids with minimal changes*.

We used Reviewed (Swiss-Prot) from [UniProt](https://www.uniprot.org/help/downloads).

- Example use cases can be found in `workspace.py`.
- Full function documentation is available in the appendix of `Tesis_Mariano_Oca.pdf`.

> \* Modifications required for:
> - `discrepancy()` in `complexity_metrics.py` (line 84)
> - `random_seq()` in `fasta_utils.py` (line 87)


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