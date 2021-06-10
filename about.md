## About the Model

Our ML analysis involved two stages: 1) fitting a cure model to predict whether redo pyeloplasty would ever occur, 2) fitting a survival model on those patients who will experience the event or are likely to (as determined by stage 1). This two-stage approach was essential to ensure that the second-stage model had an actual survival distribution. Statistical models associated with survival analysis assume that the probability of the event approaches one with sufficient measurement time. 

The cure and survival models performed well with a leave-one-out cross validation AUROC and concordance of 0.86 and 0.78, respectively. Post-selective inference showed that larger anteroposterior diameter at the second post-op follow up, and anatomical variation in the form of concurrent anomalies were significant model features predicting negative outcomes. 

More information can be found in the reference or by contacting the authors.

## Pyeloplasty code base
The code base used for this project can be found <a href="https://github.com/goldenberg-lab/Pyeloplasty">here.</a>

This repository contains the code necessary to reproduce the results and figure found in the paper: Application of Machine Learning Algorithms to Identify Patients at Risk for Recurrent UPJO after Dismembered Pyeloplasty. The scripts were run in the following order:
1. `pyl_data_prep.R`
2. `run_model.R`

The `gen_individualized.R` script allows for the calculation of a semi-parametric survival distribution based on patient-specific information which can be modified in the script.

## Reference

<b> Application of Machine Learning Algorithms to Identify Patients at Risk for Recurrent UPJO after Dismembered Pyeloplasty. </b>  <br>
<i> Drysdale E., Khondker A., Kim JK., Erdman L., Kwong JCC., Chua M., Keefe DT., Lolas M., Dos Santos J., Rickard M., Lorenzo AJ. (in preparation) 