'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'

Script that computes the performance of correctly predicting small variants (SNPs and INDELs) using the Random Forest algorithm.
It also generates a plot comparing different performances scores for each of the models tested in the project.

USAGE: python3 performances_output.py
'''

from ML_libraries_variables import *



#    C O M P U T I N G    P R O B A B I L I T I E S   ------------------------------------------------------------



datasetALL = pd.read_csv('{}/outputdat1_realvspred_randomforestrandomgrid_SW.csv'.format(os.getcwd()), sep="\t")
datasetALL.drop(['Unnamed: 0'], axis=1, inplace=True)


print()

cm = confusion_matrix(datasetALL['REAL_LABELS'],datasetALL['PREDICTED_LABELS'])

print('Confusion matrix:  \n',cm)
print()
total=sum(sum(cm))
accuracy=(cm[0,0]+cm[1,1]+cm[2,2]+cm[3,3]+cm[4,4]+cm[5,5]+cm[6,6])/total

print('Accuracy:',accuracy)



print()

for each in ['vous','polymorphism','pathogenic' ,'benign' ,'artefact' ,'riskfactor','other']:
    color_print(' {}'.format(each.upper()), color='magenta')
    print()
    real = datasetALL[datasetALL['REAL_LABELS']==each]
    predicted = datasetALL[datasetALL['PREDICTED_LABELS']==each]
    predicted_ok = datasetALL[(datasetALL['REAL_LABELS']==each) & (datasetALL['PREDICTED_LABELS']==each)]
    predicted_bad= datasetALL[(datasetALL['REAL_LABELS']!=each) & (datasetALL['PREDICTED_LABELS']==each)]
    
    print('     Real',each.upper(),len(real))
    print('     Predicted as',each.upper(),':',len(predicted))
    print('     Well predicted',each.upper(),':',len(predicted_ok))
    print('     Wrong predicted',each.upper(),':',len(predicted_bad))
    print('     Accuracy:',(round(((len(predicted_ok)/len(real))*100),2)),'%')
    

    print()
    
    real_vous = predicted_bad[(predicted_bad['REAL_LABELS']=='vous')]
    real_polimorfisme = predicted_bad[(predicted_bad['REAL_LABELS']=='polymorphism')]
    real_patogenica = predicted_bad[(predicted_bad['REAL_LABELS']=='pathogenic')]
    real_benigne = predicted_bad[(predicted_bad['REAL_LABELS']=='benign')]
    real_Artefacte = predicted_bad[(predicted_bad['REAL_LABELS']=='artefact')]
    real_riskFactor = predicted_bad[(predicted_bad['REAL_LABELS']=='riskfactor')]
    real_other = predicted_bad[(predicted_bad['REAL_LABELS']=='other')]

    if each !='vous':
        print('     The predicted is',each,'and the real is vous',len(real_vous),'    -   ', (round(((len(real_vous)/len(predicted))*100),2)),'%')
    if each !='polymorphism':
        print('     The predicted is',each,'and the real is polymorphism',len(real_polimorfisme),'    -   ', (round(((len(real_polimorfisme)/len(predicted))*100),2)),'%')
    if each !='pathogenic':
        print('     The predicted is',each,'and the real is pathogenic',len(real_patogenica),'    -   ', (round(((len(real_patogenica)/len(predicted))*100),2)),'%')
    if each !='benign':
        print('     The predicted is',each,'and the real is benign',len(real_benigne),'    -   ', (round(((len(real_benigne)/len(predicted))*100),2)),'%')
    if each !='artefact':
        print('     The predicted is',each,'and the real is artefact',len(real_Artefacte),'    -   ', (round(((len(real_Artefacte)/len(predicted))*100),2)),'%')
    if each !='riskfactor':
        print('     The predicted is',each,'and the real is riskfactor',len(real_riskFactor),'    -   ', (round(((len(real_riskFactor)/len(predicted))*100),2)),'%')
    if each !='other':
        print('     The predicted is',each,'and the real is other',len(real_other),'    -   ', (round(((len(real_other)/len(predicted))*100),2)),'%')

    print('     The predicted is ',each,'and the real is',each,len(predicted_ok),'    -   ', (round(((len(predicted_ok)/len(predicted))*100),2)),'%')
    print()


print()




#    P L O T I N G      P E R F O R M A N C E S   ----------------------------------------------------------------


fig, ax = plt.subplots(figsize=(12,12))

# A little data preparation
models = ['DTD', 'DTR', 'DTG', 'RFD', 'RFR', 'RFG', 'RFGW', 'ANND', 'ANNR', 'ANNG']
x = np.arange(len(models))*1.05

precision_benign = [0.98,0.97,0.97,0.97,0.98,0.98,0.98,0.96,0.96,0.96]
recall_benign =    [0.98,0.97,0.98,0.99,0.99,0.99,0.99,0.91,0.86,0.89] 
f1score_benign =   [0.98,0.97,0.98,0.98,0.99,0.99,0.98,0.93,0.91,0.92]
 
width=0.19

ax.barh(x - width,precision_benign,width,label='precision_benign', color='royalblue',align='center')
ax.barh(x, recall_benign, width,label='recall_benign', color='darkturquoise',align='center')
ax.barh(x + width, f1score_benign, width, label='f1score_benign', color='darkviolet',align='center')

for i, v in enumerate(precision_benign):
    ax.text(1.0,i*1.05-0.29, str(v), color='royalblue', fontsize=18) 
for i, v in enumerate(recall_benign):
    ax.text(1.0,i*1.05-0.07, str(v), color='darkturquoise', fontsize=18) 
for i, v in enumerate(f1score_benign):
    ax.text(1.0,i*1.05+0.15, str(v), color='darkviolet', fontsize=18) 


# Customise some display properties
ax.set_xlabel('performances', fontsize=18)
ax.set_ylabel('models', fontsize=18)
#ax.set_title('Model performances',y=1.04)
ax.set_xlim(0.85,1.011)
ax.set_ylim(-0.6,10.3)
ax.set_yticks(x)    # This ensures we have one tick per year, otherwise we get fewer
ax.set_yticklabels(models, rotation='horizontal')
ax.tick_params(labelsize=18)
ax.legend(bbox_to_anchor=(0,1.002 ,1,0.2), loc='lower left', ncol=3,prop={"size":18}) #(0.12,1.002 ,1,0.2)
#plt.tight_layout()
plt.savefig('model_comparisons_all_reduced_withoutaccuracies.png')
plt.close()





