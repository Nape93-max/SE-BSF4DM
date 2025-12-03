# SE-BSF4DM
Sommerfeld Effect and Bound State Formation for Dark Matter

> **Please Note:** This is the official GitHub repository for the tool **SE+BSF4DM** presented in our paper "[SE+BSF4DM: A micrOMEGAs package for Sommerfeld Effect and Bound State Formation in colored Dark Sectors](https://arxiv.org/abs/2512.02155)". The name was adapted for technical compatibility with GitHub's naming conventions.

This package extends micrOMEGAs to calculate dark matter relic density including Sommerfeld enhancement and bound state formation for QCD-colored particles.

## Installation
Follow these steps to install SE+BSF4DM:

1. Prerequisites: 
micrOMEGAs 6.0 or higher.

A CalcHEP model file for your specific model

2. Install the Core Package
Copy the entire SE_BSF directory from copy_into_Packages/ to your micrOMEGAs installation:

```bash
cp -r copy_into_Packages/SE_BSF /path/to/micrOMEGAs_6.2.4/Packages/
```

3. Model-Specific Files
For each model you want to use with SE+BSF4DM, copy the two files from copy_into_MODEL_lib/ to your model's lib/ directory:

```bash
cp copy_into_MODEL_lib/* /path/to/your/MODEL/lib/
```

Then modify these files as described in our publication.

## File Structure
copy_into_MODEL_lib/ - Model-specific files (must be copied to each model's lib/ directory)

copy_into_Packages/ - Core package files (copy SE_BSF/ to micrOMEGAs Packages/)

The Examples_Tutorial/ directory contains three models that demonstrate the simple usage of the code:

F3SuR

F3W3rd

S3Muni

These examples are provided for reference and are not required to run SE+BSF4DM.

## ATTRIBUTION
If you use **SE+BSF4DM** in your work, please cite our paper with the bibtex key:
@article{Becker:2025vgq,
    author = "Becker, Mathias and Copello, Emanuele and Napetschnig, Martin",
    title = "{SE+BSF4DM - A micrOMEGAs package for Sommerfeld Effect and Bound State Formation in colored Dark Sectors}",
    eprint = "2512.02155",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "TUM-HEP-1572/25, MITP-25-077",
    month = "12",
    year = "2025"
}
