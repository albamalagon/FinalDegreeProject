# STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS 

## Final Degree Project - Alba Malagón Márquez


- "CNV part"

  - CNV_predictions.py: Script that predicts whether a copy number variant has been correctly called or not by taking advantage of the potential presence of small variants (SNV and INDEL) in a CNV.
  - CNV_validationarrays_method.py: Script similar to 'CNV_predictions.py' but with small modifications to analyse CNVs coming from arrays, not ExomeDepth. 
  - exomedepth_and_arrays.py: Script that filters CNVs from Exome Depth and arrays with enough coverage, and stores the coincident and not coincident calls between ExomeDepth and arrays.
  - validation.py: Script to validate CNV calls between ExomeDepth, arrays and the developed method (CNV_predictions.py).
  - CNV_libraries_variables.py: all libraries and variables needed for running any of the scripts of this CNV part. 


- "ML part"
  - ML_functions.py: Script that defines all funcions to be used when predicting small variants (SNPs and INDELs) using any machine learning algorithm.
  - evaluating_algorithms.py: Script that creates and implements several machine learning algorithms to select the best one predicting small variants (SNPs and INDELs).
  - SNPandINDEL_predictions.py: Script that predicts small variants (SNPs and INDELs) using machine learning; precisely, a Random Forest algorithm.
  - train_and_test_RF.py: Script that trains and tests a Random Forest model (the best one) for predicting small variants (SNPs and INDELs).
  - performances_output.py: Script that computes the performance of correctly predicting small variants (SNPs and INDELs) using the Random Forest algorithm.
It also generates a plot comparing different performances scores for each of the models tested in the project.
  - ML_libraries_variables.py: all libraries and variables needed for running any of the scripts of this ML part. 
