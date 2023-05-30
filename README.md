# Master's thesis

[![es](https://img.shields.io/badge/lang-es-yellow.svg)](https://github.com/silvtal/TFM/blob/master/README.es.md)

This repository includes scripts from Master's Thesis "Metabolic modeling of microbial communities grown in simple energy and carbon sources" (Silvia Talavera Marcos, Master's Degree in Bioinformatics and Computational Biology, Autonomous University of Madrid). The work is downloadable [here](http://hdl.handle.net/10486/695123).

The main scripts featured here were used for annotated genome assignment, annotation, generation of metabolic models and flux balance analysis at "Leveraging phylogenetic signal to unravel microbial community function and assembly rules" by Talavera-Marcos, Parras and Aguirre de Cárcer (https://doi.org/10.21203/rs.3.rs-2272005/v1)


----------------

## General outline of the proposed pipeline

![](https://github.com/silvtal/TFM/blob/master/Anexo/Esquema.png)

## Description of the scripts

**modelado.R** includes alignment with Nucmer, model creation with CarveMe, and (optional) analysis with Smetana for a given pair of nodes. It calls the rest of the scripts, but has strict parameters hard-coded.

**annotate.R** is used for functional annotation using eggNOG-mapper and, optionally, creating a consensus model. It can also call Nucmer to start the process from the beginning. **create_compatible_database.R** is an additional script that can generate a database compatible with the version of Diamond automatically used by eggNOG-mapper from _annotate.R_.

Functions called in these scripts are located in **utils.R**.

**RefrFBA.py** is a script that performs flux balance analysis (FBA) for all models of a given phylogenetic core group (PCG). It returns the growth rates as _.csv_ files and also through standard output.

**parser.R** generates plain text reports summarizing the results of Smetana.

Finally, two Python scripts have been developed to create the consensus annotations and models: **consenso.py** creates an SBML model from other SBML models, and **consenso_EGG.py** creates a consensus annotation table from multiple annotation files.

## Input

The description of the input file format is outdated in the Appendix: now a phylogenetic tree is not required, and the file containing information about the nodes of interest is different (includes a list of the leaves for each node). Example files are included in this repository, along with a test script ("test").

## Troubleshooting

- If CarveMe/Python3.7 is installed in a Conda environment, you can activate this environment from the script itself. For example, the `carve()` function command can be modified to: `conda activate carveme && carve ",db_protein_folder`(...).

- If Smetana and/or CarveMe do not detect CPLEX or the correct version of CPLEX, the global variable `PYTHONPATH` needs to be [properly defined](https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-setting-up-python-api).
