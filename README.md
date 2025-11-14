# Primer Evaluation: A simple tool to evaluate primer properites. Based on primer3-py

This script was written to mimic the parameter output primer characteristics of NCBI Primer BLAST. However, running Primer BLAST once takes a long time, so this script was written.

## Usage (conda)

### Create conda environment

```shell
$ conda create -n primer-eval -c bioconda primer3-py
```

### Activate conda environment

```shell
$ conda activate primer-eval
```

### Run cli

```shell
$ python ./primer-eval.py <forward_primer> <reverse_primer>
```

e.g.

```$shell
$ python ./primer-eval.py "GAGTATGAGGCTTACCAGGATGGT" "GGGGTTGCGATTTTCCAGAACAA"
```