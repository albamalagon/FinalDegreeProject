'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'

Script that defines all funcions to be used when predicting small variants (SNPs and INDELs) using any machine learning algorithm.

USAGE: python3 ML_functions.py
'''


from ML_libraries_variables import *



def loading_data(alldataOK):
    '''
    function that generates a single dataset containing all SNPs and INDELs variants from the Ok.Debayan files
    '''
    color_print("### GENERATING THE WHOLE DATASET ###", color='red')
    y=1
    # creating an empty dataset where all Ok.Debayan files will be added
    column_names = ["Sample","Chr","Start","End","Ref","Alt","State","AB","DP","Filter","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","qCarrierAssay","qCarFreqInd","RunCount","qCategoria","Comment","IGV","Comment2","Procedencia","qGenFreq_536","ESP6500siv2_ALL","FreqGnomAD","NHomHemGnomAD","dbSNP","clinvar_20170905","genomicSuperDups","CV20160104","CV20160104class","X1000G_ALL","qCarVar","ExAcID","ExAcFreqObs","ExAcNRef.ExAcNObs","RunFreq","OMIM","qExFreqInd","qExVar","syn_z","mis_z","pLI","id","Num.qMales","qMales","Num.qFemales","genderBAM","qLims"]
    dataset = pd.DataFrame(columns=column_names)
    for path in alldataOK:
        color_print('Adding data {}     ({}/{})'.format(path.rsplit('/')[-1].rsplit('.')[0], y, len(alldataOK)), color='yellow')
        dataset = dataset.append(read_ods(str(path), 1, columns=column_names),ignore_index=True, sort=False)    #raises an error if the file is open!! BE CAREFUL!!
        y+=1
    return dataset


def loading_data_new(alldataOK): 
    '''
    function to read the SNPs and INDELs that need to be predicted (new data)
    '''
    color_print("### GENERATING THE WHOLE NEW DATASET ###", color='red')
    dataset = read_csv(alldataOK[0], sep="\t")
    return dataset


def labeling_and_addingpredictions(dataset, alldataGF):
    '''
    function that adds statistical predictions of the GF into the whole dataset
    '''
    newcol = ["SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "RadialSVM_pred", "LR_pred"]
    for col in newcol:
        dataset[col] = ['.'] * len(dataset) 

    i=1
    for files in alldataGF:
        genomefinal = pd.read_csv(files, sep="\t")
        identifier_GF = files.rsplit('/')[-1].rsplit('.')[0]
        color_print('Scanning GF file {}    ({}/{})'.format(identifier_GF,i,len(alldataGF)), color='cyan')
        small_dataset= dataset[dataset['Sample'] == identifier_GF]
        for index, variant in small_dataset.iterrows():
            newv=genomefinal[(variant['Chr'][-1] == genomefinal['Chr']) & (variant['Start'] == genomefinal['Start'])]
            if len(newv) > 0:
                for col in newcol:
                    dataset[col][index] = newv[col].iloc[0] 
        i+=1

    return dataset





def data_preprocessing(dataset):
    '''
    function that reduces feature dimensionality by removing uninformative features and making general simplifications of the data
    '''
    #Removing uninformative features
    dataset.drop(["Sample","Comment","qLims","qCarrierAssay","qCategoria","OMIM","IGV","Chr","Start","End","ExAcNRef.ExAcNObs","AAChange.refGene","dbSNP","CV20160104","CV20160104class","qCarVar","ExAcID",'id','qMales','qExVar','Gene.refGene','genderBAM'], axis=1, inplace=True)

    # impute missing data
    X_column_names = dataset.columns
    imputed_dataset = SimpleImputer(strategy="constant", fill_value="missing").fit_transform(dataset) 
    dataset = pd.DataFrame(imputed_dataset)
    dataset.columns = X_column_names

    # making general simplifications of the data 
    dataset.columns = (dataset.columns.str.strip().str.replace('.', '_'))
    dataset = dataset.replace(['.','na','nan','NaN','NoFet','NotFound'],'missing')
    #Procedencia
    dataset=dataset[~dataset.Procedencia.str.startswith("EDepth")]
    dataset['Procedencia'][dataset.Procedencia.str.startswith('Variant en SuperDups')] = 'missing'
    dataset=dataset[~dataset.Procedencia.str.contains("missing")] 
    #Ref and Alt
    dataset['Ref'] = dataset.Ref.replace('-','-0.0')
    dataset['Alt'] = dataset.Alt.replace('-','-0.0')
    dataset=dataset[~dataset.Alt.str.contains('0.0')] 
    #Filter
    dataset['Filter'][dataset.Filter.str.startswith('END')] = 'End'
    dataset['Filter'][dataset.Filter.str.contains('hard2validate')] = 'hard2validate'
    dataset['Filter'][dataset.Filter.str.contains('LowQual')] = 'LowQual' 
    dataset['Filter'][dataset.Filter.str.contains('LowCoverage')] = 'LowQual'
    dataset['Filter'][dataset.Filter.str.contains('SnpCluster')] = 'SnpCluster'
    dataset['Filter'][dataset.Filter.str.contains('AB_0.2')] = 'AB_0.2'
    #Clinvar
    dataset['clinvar_20170905'][dataset.clinvar_20170905.str.contains('Uncertain significance')] = 'Uncertain significance'
    dataset['clinvar_20170905'][dataset.clinvar_20170905.str.contains('enign') & dataset.clinvar_20170905.str.contains('athogenic')] = 'Uncertain significance'
    dataset['clinvar_20170905'][dataset.clinvar_20170905.str.contains('Likely benign')] = 'Likely benign'
    dataset['clinvar_20170905'][dataset.clinvar_20170905.str.contains('Benign')] = 'Benign'
    dataset['clinvar_20170905'][dataset.clinvar_20170905.str.contains('Likely pathogenic')] = 'Likely pathogenic'
    dataset['clinvar_20170905'][dataset.clinvar_20170905.str.contains('Pathogenic')] = 'Pathogenic'
    dataset['clinvar_20170905'][~dataset.clinvar_20170905.str.contains('athogenic|enign|significance')] = 'missing'
    dataset['clinvar_20170905'][dataset.clinvar_20170905.str.contains('drug response')] = 'drug response'
    #Segmental Duplications
    dataset['genomicSuperDups'][dataset.genomicSuperDups.str.contains('Score')] = 'yes'
    dataset['genomicSuperDups'] = dataset.genomicSuperDups.replace('.','no')
    #Comment2
    dataset['Comment2'][~dataset.Comment2.str.contains('Ratio_')] = 'missing' 
    dataset['Comment2'][dataset.Comment2.str.contains('Ratio_')] = dataset.Comment2.str[6:10]

    #Errors
    dataset = dataset.loc[dataset['genomicSuperDups'] != 'Check if hpalotype'] #exportar
    dataset['FreqGnomAD'] = dataset.FreqGnomAD.replace('t','missing')   
    dataset['qCarFreqInd'] = dataset.qCarFreqInd.replace('offTarget','missing') 
    dataset=dataset[~dataset.ExonicFunc_refGene.str.contains("Startloss")] 
    dataset = dataset.loc[dataset['genomicSuperDups'] != 'Check if hpalotype']

    #adding a new column that contains the length of the indel
    dataset['indel_length'] = [0] * len(dataset) 
    for indel in dataset['Ref'][dataset.Ref.str.contains('-0.0') | dataset.Alt.str.contains('-0.0')].index: 
        dataset['indel_length'][indel] = str(len(dataset['Ref'][indel])+len(dataset['Alt'][indel])-4)
        dataset['Ref'][indel] = 'indel'
        dataset['Alt'][indel] = 'indel'
    for other in dataset['Ref'][(dataset.Ref.astype(str).str.len()!=1) & (~dataset.Ref.str.contains('indel')) & (~dataset.Ref.str.contains('hba')) & (~dataset.Ref.str.contains('e7del')) & (~dataset.Alt.str.contains('TG1'))].index:
        dataset['Ref'][other] = 'other'
        dataset['Alt'][other] = 'other'

    return dataset



def minimizing_features(dataset, di):
    '''
    In case there is an unexpected category, we must decide how to handle it:
        - removing the whole variant
        - reassign the category to another one 
    '''
    for col in di.keys(): 
        dif = dataset[~dataset[col].isin(di[col])]
        if len(dif) > 0:
            print(dif)
            for each in dif[col].unique():
                color_print('There is a feature named: "{}" in the column "{}" of the dataset'.format(each, col),color='yellow')
                print('An example of a variant containing this feature looks like this:')
                print(dif.head(1))
                print('Write the correct category between these: ',di[col], 'or write "r" to remove the whole variant.')
                i = input()
                while (i != 'r') and (i not in di[col]):
                    print('The category writen (', i, ') is not between the options. Rewrite it again:')  
                    i = input()
                else:
                    if i == 'r':
                        dataset = dataset.loc[dataset[col] != each]
                    elif i in di[col]:
                        dataset[col] = dataset[col].replace(each,i)  

    return dataset


def codification(dataset, features_to_encode):

    #LABEL is a label to predict variants in y; we use the drop() function to take all other data in x. Then, we split the data.
    y = dataset.LABEL
    X = dataset.drop('LABEL',axis=1)
    dfX_cat=X[features_to_encode] #dataset containing only the categorical features
    dfX_num=X.drop(features_to_encode,axis=1) #dataset containing only the numeric features

    #dealing with missing values in numerical features
    dfX_num['DP'][dfX_num.DP.str.contains('missing')] = '-9999' #will be changed further on (by the mean)
    dfX_num['AB'][dfX_num.AB.str.contains('missing')] = '-1'

    dfX_num = dfX_num.replace('check','missing')
    dfX_num = dfX_num.replace('qCarrierPlus','missing')
    dfX_num = dfX_num.replace('CHECK','missing')
    dfX_num['syn_z'][dfX_num.syn_z.str.contains(',')] = 'missing'
    dfX_num['mis_z'][dfX_num.mis_z.str.contains(',')] = 'missing'
    dfX_num['syn_z'][dfX_num.syn_z.str.contains(';')] = 'missing'
    dfX_num['mis_z'][dfX_num.mis_z.str.contains(';')] = 'missing'
    dfX_num['pLI'][dfX_num.pLI.str.contains(',')] = 'missing'
    dfX_num['pLI'][dfX_num.pLI.str.contains(';')] = 'missing'
    dfX_num = dfX_num.replace('missing','0')

    #encoding the dataset
    #using get_dummies for the categorical features
    datasetcoded = pd.concat((dfX_num,
          pd.get_dummies(dfX_cat, dummy_na=False)),
          axis=1)

    y = y.replace('artefact','0')
    y = y.replace('benign','1')
    y = y.replace('pathogenic','2')
    y = y.replace('polymorphism','3')
    y = y.replace('riskfactor','4')
    y = y.replace('vous','5')
    y = y.replace('other','6')

    return dataset, datasetcoded, y





def compatible_dataframes(datasetcoded, datasetcodede):
    '''
    making the coded datasets compatible between them in order to make predictions
    in case some features appear in one dataset and not in the other, those features are added with '0' values on the other one
    '''

    columns_training = datasetcoded.columns
    columns_external = datasetcodede.columns

    # considering the case where the number of columns of the new dataset is lower than in the training dataset
    index=0
    for column in columns_training:
        if column not in columns_external:
            color_print('\nThe column "{}" does not appear in the new dataset'.format(column), color='yellow')
            datasetcodede.insert(index, column, [0] * len(datasetcodede))
        index+=1

    columns_training = datasetcoded.columns
    columns_external = datasetcodede.columns

    # considering the case where the number of columns of the new dataset is greater than in the training dataset
    index2=0
    for column2 in columns_external:
        if column2 not in columns_training:
            color_print('\nThe column "{}" does not appear in the training dataset'.format(column2), color='yellow')
            datasetcoded.insert(index2, column2, [0] * len(datasetcoded))
        index2+=1

    return datasetcoded, datasetcodede




def data_split(dataset, datasetcoded, y):
    '''
    splitting the X and y datasets in two: train and test
    '''

    X_train,X_test,y_train,y_test = train_test_split(datasetcoded,y,test_size=test_size, random_state=1)

    return dataset,X_train,X_test,y_train,y_test,datasetcoded





def summarizing(dataset,X_train,X_test,y_train,y_test,datasetcoded,col):
    '''
    printing some summaries of the data
    '''

    # Breakdown of the data by the class variable. The number of instances (rows) that belong to each class.
    color_print('\nNumber of elements per class:', color=col)
    print(dataset.groupby('LABEL').size())

    # Dataset dimensions
    color_print('\nDataset dimensions:', color=col)
    print(f"{dataset.shape[0]} Rows and {dataset.shape[1]} Columns")

    #first three rows of the resulting dataset
    color_print('\nHead:', color=col)
    print(dataset.head(3))

    # summarize the transformed data
    color_print('\nDimensions X and Y data:', color=col)
    color_print('X train: {}'.format(X_train.shape), color='white')
    color_print('y train: {}'.format(y_train.shape), color='white')
    color_print('X test: {}'.format(X_test.shape), color='white')
    color_print('y test: {}'.format(y_test.shape), color='white')

    #first ten rows of the resulting dataset codified
    color_print('\nHead dataset already CODED:', color=col)
    print(datasetcoded.head(10))



def summarizing_external(dataset,datasetcoded,col):
    '''
    printing some summaries of the external data
    '''
    # Breakdown of the data by the class variable. The number of instances (rows) that belong to each class.
    color_print('\nNumber of elements per class:', color=col)
    print(dataset.groupby('LABEL').size())

    # Dataset dimensions
    color_print('\nDataset dimensions:', color=col)
    print(f"{dataset.shape[0]} Rows and {dataset.shape[1]} Columns")

    #first three rows of the resulting dataset
    color_print('\nHead:', color=col)
    print(dataset.head(3))

    #first ten rows of the resulting dataset codified
    color_print('\nHead dataset already CODED:', color=col)
    print(datasetcoded.head(10))



def summarizing_new(dataset,datasetcoded,col):
    '''
    printing some summaries of the new data
    '''
    # Dataset dimensions
    color_print('\nDataset dimensions:', color=col)
    print(f"{dataset.shape[0]} Rows and {dataset.shape[1]} Columns")

    #first three rows of the resulting dataset
    color_print('\nHead:', color=col)
    print(dataset.head(3))

    #first ten rows of the resulting dataset codified
    color_print('\nHead dataset already CODED:', color=col)
    print(datasetcoded.head(10))



def evaluating_algorithm_training(algorithm, classifier, X_train, y_train, X_test, y_test, name,col):

    color_print('{}'.format(algorithm), color=col)
    if (name != 'ANN') and (name != 'ANN_r') and (name != 'ANN_g'): #considering Decision Tree and Random Forest
        start_algo_train=timeit.default_timer()
        classifier_fitted = classifier.fit(X_train,y_train)
        stop_algo_train=timeit.default_timer()
        start_algo_test=timeit.default_timer()
        #Predict the response for test dataset
        y_pred = classifier_fitted.predict(X_test)
        #Predict the response for train dataset
        y_pred_train = classifier_fitted.predict(X_train)
        stop_algo_test=timeit.default_timer()

        color_print('Time invested in training the algorithm:  {}'.format(stop_algo_train-start_algo_train), color=col)
        color_print('Time invested in testing the algorithm:  {}'.format(stop_algo_test-start_algo_test), color=col)
    
    else: #considering Artificial Neural Networks
        start_algo_test=timeit.default_timer()
        #Predict the response for test dataset
        y_pred = classifier.predict_classes(X_test)
        #Predict the response for train dataset
        y_pred_train = classifier.predict_classes(X_train)
        stop_algo_test=timeit.default_timer()
        color_print('Time invested in testing the algorithm:  {}'.format(stop_algo_test-start_algo_test), color=col)
        classifier_fitted = classifier

    
    # Model Accuracy, how often is the classifier correct?
    y_test_size = y_test.size
    y_train_size = y_train.size
    accu_train = np.sum(y_pred_train.astype(str) == y_train.astype(str))/y_train_size
    accu_test = np.sum(y_pred.astype(str) == y_test.astype(str))/y_test_size


    target_names_REAL = []
    if '0' in y_test.values:
        target_names_REAL.append('artefact')
    if '1' in y_test.values:
        target_names_REAL.append('benign')
    if '2' in y_test.values:
        target_names_REAL.append('pathogenic')        
    if '3' in y_test.values:
        target_names_REAL.append('polymorphism')
    if '4' in y_test.values:
        target_names_REAL.append('riskfactor')
    if '5' in y_test.values:
        target_names_REAL.append('vous')
    if '6' in y_test.values:
        target_names_REAL.append('other')


    print('Accuracy on train: {} %'.format(round(accu_train*100,2)))
    print('Accuracy on test: {} %'.format(round(accu_test*100,2)))

    # Printing the classification report
    color_print('\nClassification Report\n', color=col)
    print(classification_report(y_test.astype(str), y_pred.astype(str), target_names=target_names_REAL))

    classes=target_names_REAL

    # Get the confusion matrix
    cm = confusion_matrix(y_test.astype(str), y_pred.astype(str)) #, labels=classifier.classes_

    # We will store the results in a dictionary for easy access later
    per_class_accuracies = {}

    # Calculate the accuracy for each one of our classes
    for idx, cls in enumerate(classes):
        # True negatives are all the samples that are not our current GT class (not the current row) 
        # and were not predicted as the current class (not the current column)
        true_negatives = np.sum(np.delete(np.delete(cm, idx, axis=0), idx, axis=1))
        # True positives are all the samples of our current GT class that were predicted as such
        true_positives = cm[idx, idx]
        # The accuracy for the current class is ratio between correct predictions to all predictions
        per_class_accuracies[cls] = (true_positives + true_negatives) / np.sum(cm)


    color_print('\nConfusion Matrix\n', color=col)
    print(cm)
    print()
    print(per_class_accuracies)
    print()


    if (name != 'ANN') and (name != 'ANN_r') and (name != 'ANN_g'): #considering Decision Tree and Random Forest
        plot_confusion_matrix(classifier, X_test, y_test, cmap=plt.cm.Blues, display_labels=target_names_REAL,xticks_rotation=45)
        plt.tight_layout()
        plt.savefig('plot_confusionmatrix_{}.png'.format(name))
        plt.close()

        plot_confusion_matrix(classifier, X_test, y_test, cmap=plt.cm.Blues, display_labels=target_names_REAL,normalize='true',xticks_rotation=45)
        plt.tight_layout()
        plt.savefig('plot_confusionmatrix_{}_normalized.png'.format(name))
        plt.close()

        
    else: #considering Artificial Neural Networks
        
        disp=ConfusionMatrixDisplay(confusion_matrix=cm,
                                display_labels=target_names_REAL)
        disp.plot(cmap=plt.cm.Blues,xticks_rotation=45)
        plt.tight_layout()
        plt.savefig('plot_confusionmatrix_{}.png'.format(name))
        plt.close()
        

    
    # connecting predictions with real outputs
    output = X_test.copy()
    color_print('\nPREDICTIONS: ', color=col, end='')
    print(len(output))
    output['REAL_LABELS'] = y_test.astype(str)
    output['PREDICTED_LABELS'] = y_pred.astype(str)
    output['REAL_LABELS'] = output.REAL_LABELS.replace('0','artefact')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('1','benign')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('2','pathogenic')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('3','polymorphism')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('4','riskfactor')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('5','vous')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('6','other')

    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('0','artefact')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('1','benign')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('2','pathogenic')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('3','polymorphism')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('4','riskfactor')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('5','vous')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('6','other')

    print(output.head())


    outputWRONG = output.copy()
    outputWRONG = outputWRONG[~(outputWRONG['REAL_LABELS']==outputWRONG['PREDICTED_LABELS'])]
    outputWRONG.reset_index(inplace=True)

    color_print('\nWRONG PREDICTIONS: ', color=col, end='')
    print(len(outputWRONG))
    print(outputWRONG.head())

    return output, outputWRONG, classifier_fitted



def evaluating_algorithm_external(classifer_fitted, algorithm, classifier, datasetcodede, y, name, col):

    color_print('{}'.format(algorithm), color=col)

    if (name != 'ANN_external'): #considering Decision Tree and Random Forest
        #Predict the response for the whole dataset
        start_algo_test=timeit.default_timer()
        y_pred = classifer_fitted.predict(datasetcodede)
        stop_algo_test=timeit.default_timer()

        color_print('Time invested in testing the algorithm:  {}'.format(stop_algo_test-start_algo_test), color=col)
    
    else: #considering Artificial Neural Networks
        #Predict the response for the whole dataset
        y_pred = classifier.predict_classes(datasetcodede)



    target_names_REAL = []
    if '0' in y.values:
        target_names_REAL.append('artefact')
    if '1' in y.values:
        target_names_REAL.append('benign')
    if '2' in y.values:
        target_names_REAL.append('pathogenic')        
    if '3' in y.values:
        target_names_REAL.append('polymorphism')
    if '4' in y.values:
        target_names_REAL.append('riskfactor')
    if '5' in y.values:
        target_names_REAL.append('vous')
    if '6' in y.values:
        target_names_REAL.append('other')
    #if '7' in y.values:
        #target_names_REAL.append('OffTarget')


    # Printing classification report
    print('\nAccuracy: {:.2f}\n'.format(accuracy_score(y.astype(str), y_pred.astype(str))))
    color_print('\nClassification Report\n', color=col)

    print(classification_report(y.astype(str), y_pred.astype(str), target_names=target_names_REAL))

    classes=target_names_REAL

    # Get the confusion matrix
    cm = confusion_matrix(y.astype(str), y_pred.astype(str))

    # We will store the results in a dictionary for easy access later
    per_class_accuracies = {}

    # Calculate the accuracy for each one of our classes
    for idx, cls in enumerate(classes):
        # True negatives are all the samples that are not our current GT class (not the current row) 
        # and were not predicted as the current class (not the current column)
        true_negatives = np.sum(np.delete(np.delete(cm, idx, axis=0), idx, axis=1))
        # True positives are all the samples of our current GT class that were predicted as such
        true_positives = cm[idx, idx]
        # The accuracy for the current class is ratio between correct predictions to all predictions
        per_class_accuracies[cls] = (true_positives + true_negatives) / np.sum(cm)


    color_print('\nConfusion Matrix\n', color=col)
    print(cm)
    print()
    print(per_class_accuracies)

    color_print('\nMultilabel Confusion Matrix\n', color=col)
    print(multilabel_confusion_matrix(y.astype(str), y_pred.astype(str)))

    # connecting predictions with real outputs
    output = datasetcodede.copy()
    color_print('\nPREDICTIONS: ', color=col, end='')
    print(len(output))
    output['REAL_LABELS'] = y.astype(str)
    output['PREDICTED_LABELS'] = y_pred.astype(str)
    output['REAL_LABELS'] = output.REAL_LABELS.replace('0','artefact')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('1','benign')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('2','pathogenic')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('3','polymorphism')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('4','riskfactor')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('5','vous')
    output['REAL_LABELS'] = output.REAL_LABELS.replace('6','other')

    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('0','artefact')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('1','benign')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('2','pathogenic')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('3','polymorphism')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('4','riskfactor')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('5','vous')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('6','other')

    print(output.head())


    outputWRONG = output.copy()
    outputWRONG = outputWRONG[~(outputWRONG['REAL_LABELS']==outputWRONG['PREDICTED_LABELS'])]
    outputWRONG.reset_index(inplace=True)

    color_print('\nWRONG PREDICTIONS: ', color=col, end='')
    print(len(outputWRONG))
    print(outputWRONG.head())



    return output, outputWRONG



def evaluating_algorithm_new(classifer_fitted, algorithm, classifier, datasetcodede, y, name, col):

    color_print('{}'.format(algorithm), color=col)

    #Predict the response for the whole new dataset
    start_algo_test=timeit.default_timer()
    y_pred = classifer_fitted.predict(datasetcodede)
    stop_algo_test=timeit.default_timer()

    color_print('Time invested in testing the algorithm:  {}'.format(stop_algo_test-start_algo_test), color=col)
    
    output = datasetcodede.copy()
    color_print('\nPREDICTIONS: ', color=col, end='')
    print(len(output))

    output['PREDICTED_LABELS'] = y_pred.astype(str)

    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('0','artefact')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('1','benign')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('2','pathogenic')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('3','polymorphism')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('4','riskfactor')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('5','vous')
    output['PREDICTED_LABELS'] = output.PREDICTED_LABELS.replace('6','other')

    print(output.head())

    return output




def identify_zero_importance(dataset, labels, ops, n_iterations):
        '''
        Function taken from 'https://github.com/WillKoehrsen/feature-selector' and adapted to this case.
        Identify the features with zero importance according to a gradient boosting machine.
        The feature importances are averaged over `n_iterations` to reduce variance. 
        '''

        # Extract feature names
        feature_names = list(dataset.columns)
        # Convert to np array
        features = np.array(dataset)
        labels = np.array(labels).reshape((-1, ))
        # Empty array for feature importances
        feature_importance_values = np.zeros(len(feature_names))
    
        # Iterate through each fold
        for _ in range(n_iterations):
            model = RandomForestClassifier(bootstrap= False, criterion='gini', max_depth= None, max_features='auto',min_samples_leaf= 1, min_samples_split= 3, n_estimators= 1200, class_weight='balanced')
            model.fit(features, labels)

            # Record the feature importances
            feature_importance_values += model.feature_importances_ / n_iterations

        feature_importances = pd.DataFrame({'feature': feature_names, 'importance': feature_importance_values})

        # Sort features according to importance
        feature_importances = feature_importances.sort_values('importance', ascending = False).reset_index(drop = True)

        # Normalize the feature importances to add up to one
        feature_importances['normalized_importance'] = feature_importances['importance'] / feature_importances['importance'].sum()
        feature_importances['cumulative_importance'] = np.cumsum(feature_importances['normalized_importance'])

        # Extract the features with zero importance
        record_zero_importance = feature_importances[feature_importances['importance'] == 0.0]
        to_drop = list(record_zero_importance['feature'])
        feature_importances = feature_importances
        record_zero_importance = record_zero_importance
        ops['zero_importance'] = to_drop
        
        print('\n%d features with zero importance.\n' % len(ops['zero_importance']))

        # list of zero importance features
        print('The list of zero importance features:')
        print(ops['zero_importance'])
        print()

        return feature_importances, ops

def identify_low_importance(dataset, labels, cumulative_importance, feature_importances, ops):
    '''
    Function taken from 'https://github.com/WillKoehrsen/feature-selector' and adapted to this case.
    Finds the lowest importance features not needed to account for `cumulative_importance` fraction
    of the total feature importance from the random forest grid. 
    '''

    cumulative_importance = cumulative_importance
    
    # Make sure most important features are on top
    feature_importances = feature_importances.sort_values('cumulative_importance')

    # Identify the features not needed to reach the cumulative_importance
    record_low_importance = feature_importances[feature_importances['cumulative_importance'] > cumulative_importance]

    to_drop = list(record_low_importance['feature'])

    record_low_importance = record_low_importance
    ops['low_importance'] = to_drop

    print('%d features required for cumulative importance of %0.2f after one hot encoding.' % (len(feature_importances) -
                                                                        len(record_low_importance), cumulative_importance))
    print('%d features do not contribute to cumulative importance of %0.2f.\n' % (len(ops['low_importance']),
                                                                                            cumulative_importance))

    print(ops['low_importance'])

    return feature_importances, ops




def plot_feature_importances(dataset, labels, plot_n, threshold, feature_importances):
    '''
    Function taken from 'https://github.com/WillKoehrsen/feature-selector' and adapted to this case.
    Plots 'plot_n' most important features and the cumulative importance of features with a threshold.
    '''

    # Need to adjust number of features if greater than the features in the data
    if plot_n > feature_importances.shape[0]:
        plot_n = feature_importances.shape[0] - 1

    # Make a horizontal bar chart of feature importances
    plt.figure(figsize = (12, 6))
    ax = plt.subplot()

    # Need to reverse the index to plot most important on top
    features = list(reversed(list(feature_importances.index[:plot_n])))
    importances = feature_importances['normalized_importance'][:plot_n]
    plt.barh(features, importances, align = 'center', edgecolor = 'k', color='lightskyblue')

    # Set the yticks and labels
    ax.set_yticks(features)
    ax.set_yticklabels(feature_importances['feature'][:plot_n], size = 12) #, rotation=45
    ax.tick_params(labelsize=12)

    for xx,yy in zip(importances, features):
        label = "{:.5f}".format(xx)
        ax.annotate(label, # this is the text
                    (xx,yy), # this is the point to label
                    textcoords="offset points", # how to position the text
                    xytext=(30,-4), # distance from text to points (x,y)
                    ha='center', # horizontal alignment can be left, right or center
                    size=12  # rotation=90,
        )

    # Plot labeling
    plt.xlabel('Normalized Importance', size = 12)
    plt.tight_layout()
    #plt.title('Top {} Feature Importances'.format(plot_n), size = 12)
    plt.savefig('IMPORTANCES_DEFINITIVE.png')
    plt.close()

    # Cumulative importance plot
    # Index of minimum number of features needed for cumulative importance threshold
    # np.where returns the index so need to add 1 to have correct number
    plt.figure(figsize = (6, 5))
    plt.plot(list(range(1, len(feature_importances) + 1)), feature_importances['cumulative_importance'], 'b-')
    plt.xlabel('Number of Features', size = 12)
    plt.ylabel('Cumulative Importance', size = 12)
    plt.tick_params(labelsize=12)
    #plt.title('Cumulative Feature Importance', size = 12)
    importance_index = np.min(np.where(feature_importances['cumulative_importance'] > threshold))
    plt.vlines(x = importance_index + 1, ymin = 0, ymax = 1, linestyles='--', colors = 'black')
    plt.savefig('IMPORTANCES_DEFINITIVE_CUMULATIVE_threshold.png')
    plt.close()

    print('%d features required for %0.2f of cumulative importance' % (importance_index + 1, threshold))





def plot_multiclass_roc(clf, X_test, y_test, n_classes, figsize=(17, 6)):
    '''
    plotting multiclass roc curve
    '''
    y_score = clf.predict_proba(X_test)

    # structures
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # calculate dummies once
    y_test_dummies = pd.get_dummies(y_test, drop_first=False).values
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test_dummies[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # roc for each class
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('Multiclass receiver operating characteristic')
    for i in range(n_classes):
        if i == 0:
            x='artefact'
        if i == 1:
            x='benign'        
        if i == 2:
            x='pathogenic'        
        if i == 3:
            x='polymorphism'        
        if i == 4:
            x='riskfactor'        
        if i == 5:
            x='vous'       
        if i == 6:
            x='other'        

        ax.plot(fpr[i], tpr[i], label='ROC curve (area = %0.2f) for label %s' % (roc_auc[i], x))
    
    ax.legend(loc="best")
    ax.grid(alpha=.4)
    plt.savefig('ROOOOOOOC_rf.png')
    plt.close()
