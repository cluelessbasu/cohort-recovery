# CohortRecovery

This repository contains the code for the CohortRecovery algorithm, designed as a technique to produce overlapping dense subgraphs in large social networks. This is based on the paper "Finding Many Overlapping Dense Subgraphs Using Triadic Cohorts" by Sabyasachi Basu and C. Seshadhri. This work was done when the first author was a PhD student at UC Santa Cruz.

## Build

Requires GNU Make and `g++` with C++11 and OpenMP support.

```bash
make
```

## Input format

The first line contains the number of vertices and edges. Each remaining line
contains one undirected edge:

```text
<number of vertices> <number of edges>
<u> <v>
...
```

Vertex IDs must be integers from `0` to `n-1`. See `Example/small-test.txt`.

## Run

```bash
./clustering/cohorts <graph> <output-name> <epsilon> s
```

For example:

```bash
./clustering/cohorts Example/small-test.txt small-test 0.1 s
```

We recommend `epsilon = 0.1`.

The program produces the following output files in the current directory, prefaced by the graph name:

- `*-decomposition.txt`: one cluster per line
- `*-stats.txt`: cluster sizes and densities
- `*-freq.txt`: cluster-size frequencies

This implementation builds on code from
[Amazon RTRExtractor](https://github.com/amazon-science/amazon-RTRExtractor). 
See `LICENSE` and `NOTICE`.
