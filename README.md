# SE-BSF4DM
Sommerfeld Effect and Bound State Formation for Dark Matter

> **Please Note:** This is the official GitHub repository for the tool **SE+BSF4DM** presented in our papers "[Sommerfeld Effect and Bound State Formation for Dark Matter Models with Colored Mediators with SE+BSF4DM](https://arxiv.org/abs/2601.03026)" and "[Manual for SE+BSF4DM - A micrOMEGAs package for Sommerfeld Effect and Bound State Formation in colored Dark Sectors](https://arxiv.org/abs/2512.02155)". The name was adapted for technical compatibility with GitHub's naming conventions.

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
If you use **SE+BSF4DM** in your work, please cite our paper: <br>
**M. Becker, E. Copello, J. Harz and M. Napetschnig**, <br>
*"Sommerfeld Effect and Bound State Formation for Dark Matter Models with Colored Mediators with SE+BSF4DM"*, <br>
arXiv:2601.03026 [hep-ph] (2026). <br>

**M. Becker, E. Copello, J. Harz and M. Napetschnig**, <br>
*"Manual for SE+BSF4DM - A micrOMEGAs package for Sommerfeld Effect and Bound State Formation in colored Dark Sectors"*, <br>
arXiv:2512.02155 [hep-ph] (2025).
