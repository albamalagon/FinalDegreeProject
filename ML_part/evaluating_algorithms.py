'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'

Script that creates and implements several machine learning algorithms to select the best one predicting small variants (SNPs and INDELs).

USAGE: python3 evaluating_algorithms.py
'''


from ML_libraries_variables import *
from ML_functions import *

starttime = timeit.default_timer()




#    D E F I N I N G      D A T A S E T S   ----------------------------------------------------------------------------------

alldataOKtraining = OK40        # --> 40 runs
alldataGFtraining = GF40        # --> 40 runs
#alldataOKtraining = OK         # --> 21 runs
#alldataGFtraining = GF         # --> 21 runs
alldataOKexternal = OK3runs     # --> 3 runs
alldataGFexternal = GF3runs     # --> 3 runs
#alldataOKexternal = OKonerun   # --> 1 run
#alldataGFexternal = GFonerun   # --> 1 run

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

ops={}
seed = 7
numpy.random.seed(seed)







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
external_dataset = loading_data(alldataOKexternal)
external_dataset.to_csv('{}/dataset_external.csv'.format(os.getcwd()), header=True, sep="\t")
external_dataset = pd.read_csv('{}/dataset_external.csv'.format(os.getcwd()), sep="\t")
external_dataset.drop(['Unnamed: 0'], axis=1, inplace=True)


# adding statistical predictions 

#training dataset
training_dataset = labeling_and_addingpredictions(training_dataset, alldataGFtraining)
training_dataset.to_csv('{}/dataset_training_wpredictions.csv'.format(os.getcwd()), header=True, sep="\t")
training_dataset = pd.read_csv('{}/dataset_training_wpredictions.csv'.format(os.getcwd()), sep="\t")
training_dataset.drop(['Unnamed: 0'], axis=1, inplace=True)
#new dataset
external_dataset = labeling_and_addingpredictions(external_dataset, alldataGFexternal)
external_dataset.to_csv('{}/dataset_external_wpredictions.csv'.format(os.getcwd()), header=True, sep="\t")
external_dataset = pd.read_csv('{}/dataset_external_wpredictions.csv'.format(os.getcwd()), sep="\t")
external_dataset.drop(['Unnamed: 0'], axis=1, inplace=True)







#----------------------------------       STEP 2        ----------------------------------      
#                                 PREPARING THE DATASETS
#----------------------------------------------------------------------------------------- 



# LABELING THE DATA

#Adding the labels in a new column named LABEL
training_dataset['LABEL'] = training_dataset['Comment']
external_dataset['LABEL'] = external_dataset['Comment']

#joining labels with the same meaning
training_dataset['LABEL'] = training_dataset.LABEL.replace('Artefacte','artefact')
external_dataset['LABEL'] = external_dataset.LABEL.replace('Artefacte','artefact')
training_dataset['LABEL'] = training_dataset.LABEL.replace('patogenica','pathogenic')
external_dataset['LABEL'] = external_dataset.LABEL.replace('patogenica','pathogenic')
training_dataset['LABEL'] = training_dataset.LABEL.replace('polimorfisme','polymorphism')
external_dataset['LABEL'] = external_dataset.LABEL.replace('polimorfisme','polymorphism')
training_dataset['LABEL'] = training_dataset.LABEL.replace('Modificadora','riskfactor')
external_dataset['LABEL'] = external_dataset.LABEL.replace('Modificadora','riskfactor')
training_dataset['LABEL'] = training_dataset.LABEL.replace('benigne','benign')
external_dataset['LABEL'] = external_dataset.LABEL.replace('benigne','benign')
training_dataset['LABEL'] = training_dataset.LABEL.replace('Benigne','benign')
external_dataset['LABEL'] = external_dataset.LABEL.replace('Benigne','benign')
training_dataset['LABEL'] = training_dataset.LABEL.replace('sinònima','sinonima')
external_dataset['LABEL'] = external_dataset.LABEL.replace('sinònima','sinonima')



training_dataset['LABEL'] = training_dataset.LABEL.str.replace(r'(^.*ff.*$)', 'OffTarget')
external_dataset['LABEL'] = external_dataset.LABEL.str.replace(r'(^.*ff.*$)', 'OffTarget')
training_dataset['LABEL'] = np.where(~(training_dataset['LABEL'].str.contains("artefact", na=False) | (training_dataset['LABEL'].str.contains("benign", na=False)) | (training_dataset['LABEL'].str.contains("pathogenic", na=False)) | (training_dataset['LABEL'].str.contains("polymorphism", na=False)) | (training_dataset['LABEL'].str.contains("riskfactor", na=False)) | (training_dataset['LABEL'].str.contains("vous", na=False)) | (training_dataset['LABEL'].str.contains("OffTarget", na=False)) | (training_dataset['LABEL'].str.contains("sinonima", na=False))), "other", training_dataset['LABEL'])
external_dataset['LABEL'] = np.where(~(external_dataset['LABEL'].str.contains("artefact", na=False) | (external_dataset['LABEL'].str.contains("benign", na=False)) | (external_dataset['LABEL'].str.contains("pathogenic", na=False)) | (external_dataset['LABEL'].str.contains("polymorphism", na=False)) | (external_dataset['LABEL'].str.contains("riskfactor", na=False)) | (external_dataset['LABEL'].str.contains("vous", na=False)) | (external_dataset['LABEL'].str.contains("OffTarget", na=False)) | (external_dataset['LABEL'].str.contains("sinonima", na=False))), "other", external_dataset['LABEL'])

#Removing rows that contain labels we are not interested in --> OFFTARGET
# Step 1 : choose indexes to be removed
dataset_indexes_t = training_dataset[(training_dataset['LABEL'] != 'artefact') & (training_dataset['LABEL'] != 'benign') & (training_dataset['LABEL'] != 'pathogenic') & (training_dataset['LABEL'] != 'polymorphism') & (training_dataset['LABEL'] != 'vous') & (training_dataset['LABEL'] != 'riskfactor') & (training_dataset['LABEL'] != 'other')]
dataset_indexes_e = external_dataset[(external_dataset['LABEL'] != 'artefact') & (external_dataset['LABEL'] != 'benign') & (external_dataset['LABEL'] != 'pathogenic') & (external_dataset['LABEL'] != 'polymorphism') & (external_dataset['LABEL'] != 'vous') & (external_dataset['LABEL'] != 'riskfactor') & (external_dataset['LABEL'] != 'other')]
# Step 2 : remove the indexes
training_dataset = training_dataset.drop(dataset_indexes_t.index, axis=0) 
external_dataset = external_dataset.drop(dataset_indexes_e.index, axis=0) 






#----------------------------------         STEP 3        ----------------------------------      
#                                     DATA PREPROCESSING
#------------------------------------------------------------------------------------------- 


# preprocessing the data 
training_dataset = data_preprocessing(training_dataset)
external_dataset = data_preprocessing(external_dataset)


#handling new categories
training_dataset=minimizing_features(training_dataset, di)
external_dataset=minimizing_features(external_dataset, di)



training_dataset.to_csv('{}/dataset_training_definitive.csv'.format(os.getcwd()), header=True, sep="\t")
external_dataset.to_csv('{}/dataset_external_definitive.csv'.format(os.getcwd()), header=True, sep="\t")



#----------------------------------         STEP 4        ----------------------------------      
#                                     ENCODING FEATURES
#------------------------------------------------------------------------------------------- 

#categorical features
features_to_encode = ['Ref','Alt','State','Filter','Func_refGene','ExonicFunc_refGene','Procedencia','clinvar_20170905','genomicSuperDups',"SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "RadialSVM_pred", "LR_pred"]


#coding the datasets
training_dataset, datasetcoded, y = codification(training_dataset,features_to_encode)
external_dataset, datasetcodede, ye = codification(external_dataset,features_to_encode)





#----------------------------------         STEP 5        ----------------------------------      
#                 MAKING THE EXTERNAL AND TRAINING DATASETS COMPATIBLES
#------------------------------------------------------------------------------------------- 

#searching for incompatibilities between the trained dataset and the external one, and handling them
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
external_dataset['DP'][external_dataset.DP.str.contains('-9999')] = 'NaN'

X_train['DP'][X_train.DP.str.contains('NaN')] = pd.to_numeric(X_train['DP'], errors='coerce').mean()
X_test['DP'][X_test.DP.str.contains('NaN')] = pd.to_numeric(X_test['DP'], errors='coerce').mean()
external_dataset['DP'][external_dataset.DP.str.contains('NaN')] = pd.to_numeric(external_dataset['DP'], errors='coerce').mean()






#----------------------------------         STEP 8        ----------------------------------      
#                                    SUMMARIZING THE DATA
#------------------------------------------------------------------------------------------- 


color_print("### TRAINING DATASET ###", color='red')
summarizing(training_dataset,X_train,X_test,y_train,y_test,datasetcoded,col='red')

color_print("### EXTERNAL DATASET ###", color='cyan')
summarizing_external(external_dataset,datasetcodede,col='cyan')


# storing all columns appearing in the training dataset 
text_file = open("allCOLUMNStraining.txt", "w") 
text_file.write(str(datasetcoded.columns.values))
text_file.close()
# storing all columns appearing in the external dataset
text_file = open("allCOLUMNSexternal.txt", "w") 
text_file.write(str(datasetcodede.columns.values))
text_file.close()




#----------------------------------         STEP 9        ----------------------------------      
#                                  TUNNING HYPERPARAMETERS 
#------------------------------------------------------------------------------------------- 

#following the steps this article explains: https://towardsdatascience.com/hyperparameter-tuning-the-random-forest-in-python-using-scikit-learn-28d2aa77dd74

#To use RandomizedSearchCV, we first need to create a parameter grid to sample from during fitting:

# Number of trees in random forest
n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
# Number of features to consider at every split
max_features = ['auto', 'sqrt']
# Maximum number of levels in each decision tree
max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
max_depth.append(None)
# Minimum number of samples required to split a node
min_samples_split = [2, 5, 10]
# Minimum number of samples required at each leaf node
min_samples_leaf = [1, 2, 4]
# Method of selecting samples for training each tree
bootstrap = [True, False]
# The function to measure the quality of a split. Supported criteria are “gini” for the Gini impurity and “entropy” for the information gain. 
criterion = ['entropy', 'gini']

# Create the random grid   --> DECISION TREES
random_grid_DT = {"max_features": max_features,
              "max_depth": max_depth,
              'min_samples_split': min_samples_split,
              "min_samples_leaf": min_samples_leaf,
              "criterion": criterion}

# Create the random grid   --> RANDOM FOREST
random_grid_RF = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'criterion' : criterion,
               'bootstrap': bootstrap}



# Random Search to search for best hyperparameters
# Random search of parameters, using 3 fold cross validation, search across 100 different combinations, and using all available cores
#commenting this parts because it lasts a lot, and we already have know the outputs
'''
dt = DecisionTreeClassifier()
dt_random = RandomizedSearchCV(estimator = dt, param_distributions = random_grid_DT, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
dt_random.fit(X_train,y_train)
# We can view the best parameters from fitting the random search:
print('RANDOM FIT BEST PARAMETERS using DECISION TREES:')
print(dt_random.best_params_)


rf = RandomForestClassifier()
rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid_RF, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
rf_random.fit(X_train,y_train)
# We can view the best parameters from fitting the random search:
print('RANDOM FIT BEST PARAMETERS using RANDOM FOREST:')
print(rf_random.best_params_)
'''


# Grid Search with Cross Validation
# Creating the parameter grid based on the results of random search 

param_grid_dat_dt = {
    'criterion' : ['entropy'],
    'max_depth': [30, 40, 50],
    'max_features': ['auto'],
    'min_samples_split': [2, 3],
    'min_samples_leaf': [1, 2, 3]
}
param_grid_dat_rf = {
    'bootstrap': [False],
    'criterion' : ['gini'],
    'max_depth': [None],
    'max_features': ['auto'],
    'min_samples_leaf': [1, 2],
    'min_samples_split': [2, 3],
    'n_estimators': [800, 1000, 1200]
}

dt = DecisionTreeClassifier()
rf = RandomForestClassifier()
#commenting this parts because it lasts a lot, and we already have know the outputs
'''
grid_search_dt = GridSearchCV(estimator = dt, param_grid = param_grid_dat_dt, cv = 3, n_jobs = -1, verbose = 2)
grid_search_dt.fit(X_train,y_train)
print('BEST PARAMETERS GRID SEARCH using DECISION TREES:')
print(grid_search_dt.best_params_)

grid_search_rf = GridSearchCV(estimator = rf, param_grid = param_grid_dat_rf, cv = 3, n_jobs = -1, verbose = 2)
grid_search_rf.fit(X_train,y_train)
print('BEST PARAMETERS GRID SEARCH using RANDOM FOREST:')
print(grid_search_rf.best_params_)
'''

# TUNNING HYPERPARAMETERS FOR ARTIFICIAL NEURAL NETWORKS

# Function to create model, required for KerasClassifier
def create_model(optimizer, init):
	# create model
    model = Sequential()
    model.add(Dense(86, input_dim=86, activation='relu',kernel_initializer=init)) 
    model.add(Dense(70, activation='relu')) #70
    model.add(Dense(50, activation='relu')) #50
    model.add(Dense(7, activation='softmax'))
    # Compile model
    model.compile(loss='sparse_categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])   
    return model


# RANDOM SEARCH
#commenting this parts because it lasts a lot, and we already have know the outputs
'''
# grid search epochs, batch size, optimizer and init
optimizers = ['rmsprop','adam','SGD']
init = ['random_uniform','he_uniform']
epochs = [100, 150, 200, 300, 400]
batches = [20,40,55]
param_grid_ann_r = dict(optimizer=optimizers, epochs=epochs, batch_size=batches, init=init)
random_search_ann= RandomizedSearchCV(estimator=model, param_distributions=param_grid_ann_r,n_jobs=-1,cv=3, verbose = 2, random_state=42)
random_search_ann.fit(X_train, y_train)
# summarize results
print('Random Best score',random_search_ann.best_score_)
print('Random Best params',random_search_ann.best_params_)
print('Random execution time',random_search_ann.refit_time_)
'''


# GRID SEARCH
#commenting this parts because it lasts a lot, and we already have know the outputs
'''
optimizers = ['adam']
init = ['he_uniform']
epochs = [200, 300, 400]
batches = [30,35,40,45,50]
param_grid_ann_g = dict(optimizer=optimizers, epochs=epochs, batch_size=batches, init=init)
grid_search_ann = GridSearchCV(estimator = model, param_grid = param_grid_ann_g, cv = 3, n_jobs = -1, verbose = 2)
grid_search_ann.fit(X_train,y_train)
print('BEST PARAMETERS GRID SEARCH using ANN:')
print(grid_search_ann.best_params_)
'''



#TUNNING PARAMETERS OUTPUTS:

# RANDOM SEARCH

    # - DECISION TREE:
        #{'min_samples_leaf': 1, 'max_features': 'auto', 'max_depth': 40, 'criterion': 'entropy'}
    
    # - RANDOM FOREST:
        #{'n_estimators': 1000, 'min_samples_split': 2, 'min_samples_leaf': 1, 'max_features': 'auto', 'max_depth': None, 'criterion': 'gini', 'bootstrap': False}

    # - ARTIFICIAL NEURAL NETWORK:
        #{'optimizer': 'adam', 'init': 'he_uniform', 'epochs': 200, 'batch_size': 40}


# GRID SEARCH

    # - DECISION TREE:
        #{'criterion': 'entropy', 'max_depth': 40, 'max_features': 'auto', 'min_samples_leaf': 1}

    # - RANDOM FOREST:
        #{'bootstrap': False, 'criterion': 'gini', 'max_depth': None, 'max_features': 'auto', 'min_samples_leaf': 1, 'min_samples_split': 3, 'n_estimators': 1200}

    # - ARTIFICIAL NEURAL NETWORK:
        #{'batch_size': 35, 'epochs': 200, 'init': 'he_uniform', 'optimizer': 'adam'}









#----------------------------------         STEP 10        ----------------------------------      
#                                   EVALUATING ALGORITHMS
#-------------------------------------------------------------------------------------------- 

#the parameter tunning has been already performed, so we can go directly to evaluate the algorithms

#decision tree models
classifiers_dt1 = [DecisionTreeClassifier(),DecisionTreeClassifier(min_samples_leaf=1,min_samples_split= 2, max_features='auto', max_depth=40, criterion='entropy'),DecisionTreeClassifier(min_samples_leaf= 1,min_samples_split= 2, max_depth= 40, criterion= 'entropy', max_features= 'auto')]
#random forest models
classifiers_rf1 = [RandomForestClassifier(),RandomForestClassifier(n_estimators = 1000, min_samples_split= 2, min_samples_leaf= 1, max_features= 'auto', max_depth= None, criterion= 'gini', bootstrap= False), RandomForestClassifier(bootstrap= False, criterion='gini', max_depth= None, max_features='auto',min_samples_leaf= 1, min_samples_split= 3, n_estimators= 1200),RandomForestClassifier(bootstrap= False, criterion='gini', max_depth= None, max_features='auto',min_samples_leaf= 1, min_samples_split= 3, n_estimators= 1200, class_weight='balanced')]

#creating ANN models 
model1=create_model(optimizer='adam', init='he_uniform')
start_algo_train=timeit.default_timer()
model1.fit(X_train, y_train, validation_data=(X_test, y_test), epochs=400, batch_size=55)
stop_algo_train=timeit.default_timer()
color_print('Time invested in training the algorithm:  {}'.format(stop_algo_train-start_algo_train), color='red')

model2=create_model(optimizer='adam', init='he_uniform')
start_algo_train=timeit.default_timer()
model2.fit(X_train, y_train, validation_data=(X_test, y_test), epochs=400, batch_size=40)
stop_algo_train=timeit.default_timer()
color_print('Time invested in training the algorithm:  {}'.format(stop_algo_train-start_algo_train), color='red')

model3=create_model(optimizer='adam', init='he_uniform')
start_algo_train=timeit.default_timer()
model3.fit(X_train, y_train, validation_data=(X_test, y_test), epochs=400, batch_size=35)
stop_algo_train=timeit.default_timer()
color_print('Time invested in training the algorithm:  {}'.format(stop_algo_train-start_algo_train), color='red')

classifiers_ann1 = [model2, model2, model3]  

#evaluating decision tree algorithms

for i in range(len(classifiers_dt1)):
    if i == 0:
        classifier=classifiers_dt1[i]
        algorithm = '\n####        DATASET    D E C I S I O N    T R E E S        ####'
        name='decisiontree'

        outputdat1, outputWRONGdat1, classifer_fitted_dt1 = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test,  name, col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")

    elif i == 1:
        classifier=classifiers_dt1[i]
        algorithm = '\n####        DATASET    D E C I S I O N    T R E E S    R A N D O M        ####'
        name='decisiontreerandom'

        outputdat1, outputWRONGdat1, classifer_fitted_dt1_random = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")

    elif i == 2:
        classifier=classifiers_dt1[i]
        algorithm = '\n####        DATASET    D E C I S I O N    T R E E S    R A N D O M      G R I D        ####'
        name='decisiontreerandomgrid'

        outputdat1, outputWRONGdat1, classifer_fitted_dt1_random_grid = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")

#evaluating random forest algorithms
#adding also some external evaluations after we have observe Random Forest is the best performing algorithm

for i in range(len(classifiers_rf1)):
    
    if i == 0:
        classifier=classifiers_rf1[i]
        algorithm = '\n####        DATASET    R A N D O M     F O R E S T        ####'
        name='randomforest'

        outputdat1, outputWRONGdat1, classifer_fitted_rf1 = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name,col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
    
    if i == 1:
        classifier=classifiers_rf1[i]
        algorithm = '\n####        DATASET    R A N D O M     F O R E S T    R A N D O M        ####'
        name='randomforestrandom'

        outputdat1, outputWRONGdat1, classifer_fitted_rf1_random = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
    
        algorithm = '\n####        DATASET  EXTERNAL    R A N D O M     F O R E S T    R A N D O M        ####'
        name='RFRANDOM'
        outputdat1e, outputWRONGdat1e = evaluating_algorithm_external(classifer_fitted_rf1_random, algorithm, classifier, datasetcodede, ye, name, col='cyan')
        outputdat1e.to_csv('{}/outputdat1e_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1e.to_csv('{}/outputdat1e_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
    

    if i == 2:
        classifier=classifiers_rf1[i]
        algorithm = '\n####        DATASET    R A N D O M     F O R E S T    R A N D O M      G R I D        ####'
        name='randomforestrandomgrid'

        outputdat1, outputWRONGdat1, classifer_fitted_rf1_random_grid = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
        
        algorithm = '\n####        DATASET  EXTERNAL    R A N D O M     F O R E S T    R A N D O M     G R I D    ####'
        name='RFRANDOMGRID'
        outputdat1e, outputWRONGdat1e = evaluating_algorithm_external(classifer_fitted_rf1_random_grid, algorithm, classifier, datasetcodede, ye, name, col='cyan')
        outputdat1e.to_csv('{}/outputdat1e_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1e.to_csv('{}/outputdat1e_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
        

    if i == 3: #adding sample weight

        classifier=classifiers_rf1[i]
        algorithm = '\n####        DATASET    R A N D O M     F O R E S T    R A N D O M      G R I D     SW    ####'
        name='randomforestrandomgrid_SW'

        outputdat1, outputWRONGdat1, classifer_fitted_rf1_random_grid_SW = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='red')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
        
        algorithm = '\n####        DATASET  EXTERNAL    R A N D O M     F O R E S T    R A N D O M     G R I D      SW  ####'
        name='RFRANDOMGRID_SW'
        outputdat1e, outputWRONGdat1e = evaluating_algorithm_external(classifer_fitted_rf1_random_grid_SW, algorithm, classifier, datasetcodede, ye, name, col='cyan')
        outputdat1e.to_csv('{}/outputdat1e_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1e.to_csv('{}/outputdat1e_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")



#evaluating Artificial Neural Networks algorithms
#adding also an external test

for i in range(len(classifiers_ann1)):
    
    if i == 0:
        classifier=classifiers_ann1[i]
        algorithm = '\n####        DATASET    A R T I F I C I A L    N E U R A L     N E T W O R K S    ####'
        name='ANN'

        outputdat1, outputWRONGdat1, classifer_fitted_ANN = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='magenta')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
        
        algorithm = '\n####        DATASET   EXTERNAL    A R T I F I C I A L    N E U R A L     N E T W O R K S    ####'
        name='ANN_external'

        outputdat1, outputWRONGdat1 = evaluating_algorithm_external(classifer_fitted_ANN,algorithm, classifier, datasetcodede, ye, name, col='cyan')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")
    
    
    if i == 1:
        classifier=classifiers_ann1[i]
        algorithm = '\n####        DATASET    A R T I F I C I A L    N E U R A L     N E T W O R K S     R A N D O M ####'
        name='ANN_r'

        outputdat1, outputWRONGdat1, classifer_fitted_ANN = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='magenta')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")

    
    if i == 2: 
        classifier=classifiers_ann1[i]
        algorithm = '\n####        DATASET    A R T I F I C I A L    N E U R A L     N E T W O R K S      G R I D   ####'
        name='ANN_g'

        outputdat1, outputWRONGdat1, classifer_fitted_ANN = evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name, col='magenta')
        outputdat1.to_csv('{}/outputdat1_realvspred_{}.csv'.format(os.getcwd(),name), header=True, sep="\t")
        outputWRONGdat1.to_csv('{}/outputdat1_realvspred_{}_WRONG.csv'.format(os.getcwd(),name), header=True, sep="\t")


     






#----------------------------------         STEP 11        ----------------------------------      
#                                         SOME PLOTS
#-------------------------------------------------------------------------------------------- 


# PLOTING FEATURE IMPORTANCES

feature_importances = pd.read_csv('{}/feature_importances.csv'.format(os.getcwd()), sep="\t")
feature_importances.drop(['Unnamed: 0'], axis=1, inplace=True)


plot_feature_importances(datasetcoded, y, num_feature_importance, threshold, feature_importances)





# PLOTING one tree from the forest
'''
features_list = datasetcoded.columns

# Pull out one tree from the forest
tree = classifer_fitted_rf1_random_grid.estimators_[5]

# DOT data
dot_data = export_graphviz(tree, out_file=None, 
                                feature_names=features_list,
                                class_names=['artefact','benign','pathogenic','polymorphism','riskfactor','vous', 'other'],
                                filled=True,
                                rounded=True)

# Draw graph
graph = graphviz.Source(dot_data, format="png") 
graph.render("randomforestgrid_graphivz")
'''


# PLOTING a reduced tree from the forest
'''
rf_small=RandomForestClassifier(bootstrap= False, criterion='gini', max_depth= 3, max_features='auto',min_samples_leaf= 1, min_samples_split= 3, n_estimators= 10)
rf_fitted_small = rf_small.fit(X_train,y_train)

tree_small = rf_fitted_small.estimators_[0]
# DOT data
dot_data_small = export_graphviz(tree_small, out_file=None, 
                                feature_names=features_list,
                                class_names=['artefact','benign','pathogenic','polymorphism','riskfactor','vous', 'other'],
                                filled=True,
                                rounded=True)

# Draw graph
graph_small = graphviz.Source(dot_data_small, format="png") 
graph_small.render("randomforestgrid_graphivz_small")
'''




stoptime = timeit.default_timer()
color_print('\nTime:', color='red')
print(stoptime - starttime)



