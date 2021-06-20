# STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS 

## Final Degree Project - Alba Malagón Márquez

_The motivation behind this project is to reduce the time spent by the experts analysing the variants through two related strategies. The first, the improvement in the detection of copy number variants by assessing their credibility. The second, improved efficiency in the tertiary analysis of NGS data by implementing machine learning tools that can predict their clinical importance. In essence, both parts are aimed at improving the efficiency of NGS processes in the interpretation part._


### FILES

- CNV part

  - [CNV_predictions.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/CNV_part/CNV_predictions.py): Script that predicts whether a copy number variant has been correctly called or not by taking advantage of the potential presence of small variants (SNV and INDEL) in a CNV.
  - [CNV_validationarrays_method.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/CNV_part/CNV_validation_method.py): Script similar to 'CNV_predictions.py' but with small modifications to analyse CNVs coming from arrays, not ExomeDepth. 
  - [exomedepth_and_arrays.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/CNV_part/exomedepth_and_arrays.py): Script that filters CNVs from Exome Depth and arrays with enough coverage, and stores the coincident and not coincident calls between ExomeDepth and arrays.
  - [validation.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/CNV_part/validation.py): Script to validate CNV calls between ExomeDepth, arrays and the developed method (CNV_predictions.py).
  - [CNV_libraries_variables.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/CNV_part/CNV_libraries_variables.py): all libraries and variables needed for running any of the scripts of this CNV part. 


- "ML part"
  - [ML_functions.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/ML_part/ML_functions.py): Script that defines all funcions to be used when predicting small variants (SNPs and INDELs) using any machine learning algorithm.
  - [evaluating_algorithms.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/ML_part/evaluating_algorithms.py): Script that creates and implements several machine learning algorithms to select the best one predicting small variants (SNPs and INDELs).
  - [SNPandINDEL_predictions.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/ML_part/SNPandINDEL_predictions.py): Script that predicts small variants (SNPs and INDELs) using machine learning; precisely, a Random Forest algorithm.
  - [train_and_test_RF.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/ML_part/train_and_test_RF.py): Script that trains and tests a Random Forest model (the best one) for predicting small variants (SNPs and INDELs).
  - [performances_output.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/ML_part/performances_output.py): Script that computes the performance of correctly predicting small variants (SNPs and INDELs) using the Random Forest algorithm.
It also generates a plot comparing different performances scores for each of the models tested in the project.
  - [ML_libraries_variables.py](https://github.com/albamalagon/FinalDegreeProject/blob/main/ML_part/ML_libraries_variables.py): all libraries and variables needed for running any of the scripts of this ML part. 


### USAGE

All scripts are executed in the same way:
```
python3 filename.py
```
```
