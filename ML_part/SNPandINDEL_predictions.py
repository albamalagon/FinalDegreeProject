'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'

Script that predicts small variants (SNPs and INDELs) using machine learning; precisely, a Random Forest algorithm.

USAGE: python3 SNPandINDEL_predictions.py
'''


from ML_libraries_variables import *
from ML_functions import *

starttime = timeit.default_timer()



color_print('Write which run number you desire to predict:',color='red')

num_run = input()


#    D E F I N I N G      D A T A S E T S   ----------------------------------------------------------------------------------

alldataOKtraining = OK40         # --> 40 runs
alldataGFtraining = GF40         # --> 40 runs
alldataOKnew= glob.glob(f"{os.getcwd()}/run{num_run}/OK/*.qCarrier.ALLsamples.AGRUPANT_TAULES.qAnotat.OK.debayan.txt")
alldataGFnew = glob.glob(f"{os.getcwd()}/run{num_run}/GF/*.genome.FINAL.txt")

#    D E F I N I N G      F E A T U R E S    ----------------------------------------------------------------------------------

#these are the only categories each feature can contain
ref=['A','C','G','T','indel','other','missing','TT','HBA2']
alt=['A','C','G','T','indel','other','TG10','TG11','TG13','TG12','missing','-0.0']
state=['hom','het','CHECK', 'unknown','missing']
filterr=['PASS','AB_0.2','LowQual','SnpCluster','End','hard2validate','missing']
funcrefgene=['exonic','splicing','intronic','ncRNA_exonic','ncRNA_intronic','ncRNA_splicing','UTR3','missing','upstream','exonic;splicing','ncRNA_exonic;splicing']
procedencia=['No1en25','ClinVar_NoClinVar','NoCVar20160104','missing','SMN1file','TGfile','HBAbedFile']
exonicfuncrefgene=['nonsynonymous_SNV','splicing','frameshift_insertion','stoploss', 
 'frameshift_deletion','nonframeshift_insertion','duplication','deletion',
 'nonframeshift_deletion','stopgain','PseudoGen_deletion',
 'PseudoGen_duplication','c.1210-34TG(10)T(5)','c.1210-34TG(11)T(5)',
 'c.1210-34TG(13)T(5)','c.1210-34TG(12)T(5)','missing','unknown',
 'synonymous_SNV','nonframeshift_substitution','frameshift_substitution','exonic;splicing', 'ncRNA_splicing','ncRNA_exonic;splicing']        
clinvar=['Uncertain significance','missing','Pathogenic','Likely pathogenic',
 'Benign','Likely benign','drug response']
genomicsuperdups=['yes', 'no', 'missing']

di = {'Ref': ref, 'Alt': alt, 'State' : state, 'Filter' : filterr, 'Func_refGene' : funcrefgene, 'ExonicFunc_refGene': exonicfuncrefgene,
      'Procedencia': procedencia, 'clinvar_20170905' : clinvar, 'genomicSuperDups': genomicsuperdups}

# O T H E R      V A R I A B L E S      --------------------------------------------------------------------------------------
seed=7
numpy.random.seed(seed)
ops={}







#----------------------------------       STEP 1       ----------------------------------      
#                                    LOADING THE DATA
#----------------------------------------------------------------------------------------- 


# generating the whole datasets

#training dataset
training_dataset = loading_data(alldataOKtraining)
training_dataset.to_csv('{}/dataset_training.csv'.format(os.getcwd()), header=True, sep="\t")
training_dataset = pd.read_csv('{}/dataset_training.csv'.format(os.getcwd()), sep="\t")
training_dataset.drop(['Unnamed: 0'], axis=1, inplace=True)
#new dataset
new_dataset = loading_data_new(alldataOKnew)
new_dataset.to_csv('{}/dataset_new_{}.csv'.format(os.getcwd(),num_run), header=True, sep="\t")
new_dataset = pd.read_csv('{}/dataset_new_{}.csv'.format(os.getcwd(),num_run), sep="\t")
new_dataset.drop(['Unnamed: 0'], axis=1, inplace=True)

# adding statistical predictions 

#training dataset
training_dataset = labeling_and_addingpredictions(training_dataset, alldataGFtraining)
training_dataset.to_csv('{}/dataset_training_wpredictions.csv'.format(os.getcwd()), header=True, sep="\t")
training_dataset = pd.read_csv('{}/dataset_training_wpredictions.csv'.format(os.getcwd()), sep="\t")
training_dataset.drop(['Unnamed: 0'], axis=1, inplace=True)
#new dataset
new_dataset = labeling_and_addingpredictions(new_dataset, alldataGFnew)
new_dataset.to_csv('{}/dataset_new_wpredictions_{}.csv'.format(os.getcwd(),num_run), header=True, sep="\t")
new_dataset = pd.read_csv('{}/dataset_new_wpredictions_{}.csv'.format(os.getcwd(),num_run), sep="\t")
new_dataset.drop(['Unnamed: 0'], axis=1, inplace=True)



whole_training_dataset = pd.read_csv('{}/dataset_training_wpredictions.csv'.format(os.getcwd()), sep="\t")
whole_new_dataset = pd.read_csv('{}/dataset_new_wpredictions_{}.csv'.format(os.getcwd(),num_run), sep="\t")





#----------------------------------       STEP 2        ----------------------------------      
#                             PREPARING THE TRAINING DATASET
#----------------------------------------------------------------------------------------- 


# LABELING THE DATA

#adding the labels in a new column named LABEL
training_dataset['LABEL'] = training_dataset['Comment']
new_dataset['LABEL'] = new_dataset['Comment']

#joining labels with the same meaning
training_dataset['LABEL'] = training_dataset.LABEL.replace('Artefacte','artefact')
training_dataset['LABEL'] = training_dataset.LABEL.replace('patogenica','pathogenic')
training_dataset['LABEL'] = training_dataset.LABEL.replace('polimorfisme','polymorphism')
training_dataset['LABEL'] = training_dataset.LABEL.replace('Modificadora','riskfactor')
training_dataset['LABEL'] = training_dataset.LABEL.replace('benigne','benign')
training_dataset['LABEL'] = training_dataset.LABEL.replace('Benigne','benign')
training_dataset['LABEL'] = training_dataset.LABEL.replace('sinònima','sinonima')

training_dataset['LABEL'] = training_dataset.LABEL.str.replace(r'(^.*ff.*$)', 'OffTarget')
training_dataset['LABEL'] = np.where(~(training_dataset['LABEL'].str.contains("artefact", na=False) | (training_dataset['LABEL'].str.contains("benign", na=False)) | (training_dataset['LABEL'].str.contains("pathogenic", na=False)) | (training_dataset['LABEL'].str.contains("polymorphism", na=False)) | (training_dataset['LABEL'].str.contains("riskfactor", na=False)) | (training_dataset['LABEL'].str.contains("vous", na=False)) | (training_dataset['LABEL'].str.contains("OffTarget", na=False)) | (training_dataset['LABEL'].str.contains("sinonima", na=False))), "other", training_dataset['LABEL'])

#Removing rows that contain labels we are not interested in --> OFFTARGET
# Step 1 : choose indexes to be removed
dataset_indexes_t = training_dataset[(training_dataset['LABEL'] != 'artefact') & (training_dataset['LABEL'] != 'benign') & (training_dataset['LABEL'] != 'pathogenic') & (training_dataset['LABEL'] != 'polymorphism') & (training_dataset['LABEL'] != 'vous') & (training_dataset['LABEL'] != 'riskfactor') & (training_dataset['LABEL'] != 'other')]
# Step 2 : remove the indexes
training_dataset = training_dataset.drop(dataset_indexes_t.index, axis=0) 







#----------------------------------         STEP 3        ----------------------------------      
#                                     DATA PREPROCESSING
#------------------------------------------------------------------------------------------- 


# preprocessing the data 
training_dataset = data_preprocessing(training_dataset)
new_dataset = data_preprocessing(new_dataset)


#handling new categories
training_dataset=minimizing_features(training_dataset, di)
new_dataset=minimizing_features(new_dataset, di)


training_dataset.to_csv('{}/dataset_training_definitive.csv'.format(os.getcwd()), header=True, sep="\t")
new_dataset.to_csv('{}/dataset_new_definitive_{}.csv'.format(os.getcwd(),num_run), header=True, sep="\t")







#----------------------------------         STEP 4        ----------------------------------      
#                                     ENCODING FEATURES
#------------------------------------------------------------------------------------------- 

#categorical features
features_to_encode = ['Ref','Alt','State','Filter','Func_refGene','ExonicFunc_refGene','Procedencia','clinvar_20170905','genomicSuperDups',"SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "RadialSVM_pred", "LR_pred"]


#coding the datasets
training_dataset, datasetcoded, y = codification(training_dataset,features_to_encode)
new_dataset, datasetcodede, ye = codification(new_dataset,features_to_encode)





#----------------------------------         STEP 5        ----------------------------------      
#                 MAKING THE EXTERNAL AND TRAINING DATASETS COMPATIBLES
#------------------------------------------------------------------------------------------- 

#searching for incompatibilities between the trained dataset and the new one, and handling them
datasetcoded, datasetcodede = compatible_dataframes(datasetcoded, datasetcodede)





#----------------------------------         STEP 6        ----------------------------------      
#                    FEATURE SELECTION: REMOVING LOW IMPORTANCES FEATURES
#------------------------------------------------------------------------------------------- 


#selecting zero and low importance features 
#in this case it has been done with the best perfomance model: random forest random with specific grid weighted
feature_importances, ops = identify_zero_importance(datasetcoded, y, ops, n_iterations)

feature_importances, ops = identify_low_importance(datasetcoded, y, cumulative_importance, feature_importances, ops)

datasetcoded=datasetcoded.drop(ops['low_importance'], axis=1)
datasetcodede=datasetcodede.drop(ops['low_importance'], axis=1)


f = open("dict.txt","w")
f.write( str(ops) )
f.close()


feature_importances.to_csv('{}/feature_importances.csv'.format(os.getcwd()), header=True, sep="\t")








#----------------------------------         STEP 7        ----------------------------------      
#                             SPLITTING THE DATA (training only)
#------------------------------------------------------------------------------------------- 


#Let’s split this data into training and test. 
training_dataset,X_train,X_test,y_train,y_test,datasetcoded = data_split(training_dataset, datasetcoded, y)



# considering missing values: case of mean

X_train['DP'][X_train.DP.str.contains('-9999')] = 'NaN'
X_test['DP'][X_test.DP.str.contains('-9999')] = 'NaN'
new_dataset['DP'][new_dataset.DP.str.contains('-9999')] = 'NaN'

X_train['DP'][X_train.DP.str.contains('NaN')] = pd.to_numeric(X_train['DP'], errors='coerce').mean()
X_test['DP'][X_test.DP.str.contains('NaN')] = pd.to_numeric(X_test['DP'], errors='coerce').mean()
new_dataset['DP'][new_dataset.DP.str.contains('NaN')] = pd.to_numeric(new_dataset['DP'], errors='coerce').mean()






#----------------------------------         STEP 8        ----------------------------------      
#                                    SUMMARIZING THE DATA
#------------------------------------------------------------------------------------------- 


color_print("### TRAINING DATASET ###", color='red')
summarizing(training_dataset,X_train,X_test,y_train,y_test,datasetcoded,col='red')

color_print("### NEW DATASET ###", color='cyan')
summarizing_new(new_dataset,datasetcodede,col='cyan')


# storing all columns appearing in the training dataset 
text_file = open("allCOLUMNStraining.txt", "w") 
text_file.write(str(datasetcoded.columns.values))
text_file.close()
# storing all columns appearing in the external dataset
text_file = open("allCOLUMNSexternal.txt", "w") 
text_file.write(str(datasetcodede.columns.values))
text_file.close()






#----------------------------------         STEP 10        ----------------------------------      
#                                   EVALUATING ALGORITHMS
#-------------------------------------------------------------------------------------------- 

#the parameter tunning has been already performed, so we can go directly to evaluate the algorithm
classifiers_rf1 = [RandomForestClassifier(bootstrap= False, criterion='gini', max_depth= None, max_features='auto',min_samples_leaf= 1, min_samples_split= 3, n_estimators= 1200)]


for i in range(len(classifiers_rf1)):

    if i == 0:
        classifier=classifiers_rf1[i]
        algorithm = '\n####        DATASET    R A N D O M     F O R E S T    R A N D O M      G R I D     W E I G H T E D     ####'
        name='randomforestrandomgridw'

        outputdat1, outputWRONGdat1, classifer_fitted_rf1_random_grid = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
        
        algorithm = '\n####        D A T A S E T     N E W    with   R A N D O M     F O R E S T    R A N D O M     G R I D     W E I G H T E D     ####'
        name='randomforestrandomgridw_new'
        outputnew = evaluating_algorithm_new(classifer_fitted_rf1_random_grid, algorithm, classifier, datasetcodede, ye, name, col='cyan')
        outputnew.to_csv('{}/output_{}.csv'.format(os.getcwd(),num_run), header=True, sep="\t")
        



# Look at parameters used by the best model
print('Parameters of the outperforming model:\n')
print(classifer_fitted_rf1_random_grid.get_params())


#joining the predicted variants with the ones not used

outputdat1 = pd.read_csv('{}/outputdat1_realvspred_randomforestrandomgridw.csv'.format(os.getcwd()), sep="\t")
outputnew = pd.read_csv('{}/output_{}.csv'.format(os.getcwd(),num_run), sep="\t")

allpredictions = whole_training_dataset.merge(outputdat1, how = 'left' , on = ['Unnamed: 0'], indicator=True)
allpredictionse = whole_new_dataset.merge(outputnew, how = 'left' , on = ['Unnamed: 0'], indicator=True)

allpredictions.columns = (allpredictions.columns.str.strip().str.replace('_x', ''))
allpredictionse.columns = (allpredictionse.columns.str.strip().str.replace('_x', ''))

delete_col=[]
delete_cole=[]

for c in allpredictions.columns:
    if (c not in whole_training_dataset.columns) and (c !='REAL_LABELS') and (c !='PREDICTED_LABELS'):
        delete_col.append(c)
    if c[-2:] == '_y':
        delete_col.append(c)

for ce in allpredictionse.columns:
    if (ce not in whole_new_dataset.columns) and (ce !='PREDICTED_LABELS'):
        delete_cole.append(ce)
    if ce[-2:] == '_y':
        delete_cole.append(ce)

allpredictions.drop(delete_col, axis=1, inplace=True)
allpredictions.drop(['Unnamed: 0'], axis=1, inplace=True)
allpredictionse.drop(delete_cole, axis=1, inplace=True)
allpredictionse.drop(['Unnamed: 0'], axis=1, inplace=True)



allpredictions.to_csv('{}/allpredictions.csv'.format(os.getcwd()), header=True, sep="\t")
allpredictionse.to_csv('{}/allpredictions_{}.csv'.format(os.getcwd(),num_run), header=True, sep="\t")



color_print('\nDataset dimensions:', color='green')
print(f"{allpredictions.shape[0]} Rows and {allpredictions.shape[1]} Columns")
print(allpredictions.head(15))
color_print('\nDataset dimensions:', color='green')
print(f"{allpredictionse.shape[0]} Rows and {allpredictionse.shape[1]} Columns")
print(allpredictionse.head(15))










#----------------------------------         STEP 11        ----------------------------------      
#                                       PLOTING RESULTS
#-------------------------------------------------------------------------------------------- 


# ploting multiclass roc curve
plot_multiclass_roc(classifer_fitted_rf1_random_grid, X_test, y_test, n_classes=7, figsize=(6, 5))

#ploting feature importances
feature_importances = pd.read_csv('{}/feature_importances.csv'.format(os.getcwd()), sep="\t")
feature_importances.drop(['Unnamed: 0'], axis=1, inplace=True)
plot_feature_importances(datasetcoded, y, num_feature_importance, threshold, feature_importances)









stoptime = timeit.default_timer()
color_print('\nTime:', color='red')
print(stoptime - starttime)



