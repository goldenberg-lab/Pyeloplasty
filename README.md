## Pyeloplasty code base

This repository contains the code necessary to reproduce the results and figure found in the paper: Application of Machine Learning Algorithms to Identify Patients at Risk for Recurrent UPJO after Dismembered Pyeloplasty. The scripts were run in the following order:

1. `pyl_data_prep.R`
2. `run_model.R`


The `gen_individualized.R` script allows for the calculation of a semi-parametric survival distribution based on patient-specific information which can be modified in the script.


