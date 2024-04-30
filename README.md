# HighFold: accurately predicting structures of cyclic peptides and complexes with head-to-tail and disulfide bridge constraints

HighFold is a novol framework for cyclic peptide structure prediction. Based on head-to-tail and disulfide bridge constraints, a Cyclic Position Offset Encoding Matrix (CycPOEM) is constructed and feed into the AlphaFold model. 

## Advantanges of HighFold

- Cyclic peptide structure prediction for both monomers and complexes.
- Multiple cyclizations such as head-to-tail and disulfide bridges.
- High effectiveness and efficiency.

## Code Hierarchy

```shell
bashCopy highfold/
├── alphafold/          # core codes from alphafold
├── colabfold/          # colabfold implement
├── utils/              # construction of CycPOEM, combination of disufide bridges and metrics
├── LICENSE             # license file
└── README.md           # readme
└── HighFold_data       # HighFold_data is the dataset used in our work.
```

## Instructions

### Installation

You can install the ColabFold by the script LocalColabFold (details on https://github.com/YoshitakaMo/localcolabfold) at first, and then copy the source codes of HighFold into the installed ColabFold preject.

### Data Preparation

The amino acid sequence should be processed in fasta format.

```
>fasta
SAKIDNLD:
SSPGIWLDCTHLEGKVILVAVHVASGYIEAVIPAETGQETAYFLLLAGRWPVKTHDNGSNFTSTTVKAACWWAGIQEDGIPYNPQSQGVIESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNHKRKGYSAGERIVIIATDIQ:
SSPGIWLDCTHLEGKVILVAVHVASGYIEAVIPAETGQETAYFLLLAGRWPVKTHDNGSNFTSTTVKAACWWAGIQEDGIPYNPQSQGVIESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNHKRKGYSAGERIVIIATDIQ
```

### Structure Prediction

To run the prediction, type

```sh
colabfold_batch --templates --amber --model-type alphafold2 file_input path_output [args]
```

## Version History

- v1.0.0 (2023-05-18):

