'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'

Script that filters CNVs from Exome Depth and arrays with enough coverage, and stores the coincident and not coincident calls between ExomeDepth and arrays.

USAGE: python3 exomedepth_and_arrays.py
'''

from CNV_libraries_variables import *

starttime = timeit.default_timer()



#    D E F I N I N G    V A R I A B L E S    ----------------------------------------------------------------------

alldataGFval = GFval 
alldataEDval = EDval 
alldataarray = array


#    R E A D I N G    F I L E S   ----------------------------------------------------------------------------------

#reading arrays designs
design60_37 = pd.read_csv(design[0], sep='\t')
design60_24 = pd.read_csv(design[1], sep='\t')
designPrePost = pd.read_csv(design[2], sep='\t')
design400 = pd.read_csv(design[3], sep='\t')
design180 = pd.read_csv(design[4], sep='\t')
design1M = pd.read_csv(design[5], sep='\t')
#reading samples with array and exoma
exomaa = pd.read_csv(exoma[0],sep=';')
#reading the pseudogenes dataset
pseudogenes = pd.read_csv(pseudogenes[0], sep='\t')
#reading the segmental duplications
segdups = pd.read_csv(SD[0], sep='\t')
#reading exoma data
capexoma = pd.read_csv(captura_ex[0],sep='\t')


color_print("### GENERATING THE WHOLE DATASET ###", color='red')
y=1
# creating an empty dataset where all Ok.Debayan files will be added
column_names = ['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'OMIM', 'OMIMFeb2017', 'NqExDels',
       'NqExDups', 'qExDels', 'qExDups','Sample']
dataset = pd.DataFrame(columns=column_names)


#Adding ExomeDepth CNVs data
for path in alldataEDval:
    color_print('Adding EXOME DEPTH data {}     ({}/{})'.format(path.rsplit('/')[-1].rsplit('_')[0], y, len(alldataEDval)), color='yellow')
    variable_length=len(dataset)
    dataset = dataset.append(pd.read_csv(str(path), sep="\t"),ignore_index=True, sort=False)    #raises and error if the file is open!! BE CAREFUL!
    dataset.loc[variable_length:len(dataset),['Sample']] = path.rsplit('/')[-1].rsplit('_')[0]
    y+=1

dataset.to_csv('{}/datasetCNV_VALIDATION.csv'.format(os.getcwd()), header=True, sep="\t")

dataset = pd.read_csv('{}/datasetCNV_VALIDATION.csv'.format(os.getcwd()), sep="\t")

dataset.drop(['OMIM', 'OMIMFeb2017', 'NqExDels',
       'NqExDups', 'qExDels', 'qExDups', 'Unnamed: 0'], axis=1, inplace=True)


#    C O V E R E D   C N V S   E X O M E   D E P T H   -------------------------------------------------------------


#creating three new columns to keep track of the probes, coverage and the array type of each CNV
dataset['PROBES'] = [0] * len(dataset) 
dataset['COVERAGE'] = 'False'
dataset["array_type"] = ['.'] * len(dataset) 

# analysing each CNV individually
for index, cnv in dataset.iterrows():
    chro = cnv['chr']
    start = cnv['start']
    end = cnv['end']
    typee = cnv['type']
    allchr='chr'+str(cnv['chr'])

    #adding the array type corresponding to the current sample
    array=exomaa[(exomaa['Mostra'] == dataset['Sample'][index])] 
    dataset['array_type'][index] = array['Tipus_prova_array'].iloc[0]


    #extracting the number of probes the CNV contains taking into account its design type
    if (dataset['array_type'][index] == 'qChip Post') or (dataset['array_type'][index] == 'qChip Pre v1.1 Complete') or (dataset['array_type'][index] == 'qChip Pre v1.1 Targeted'):
        num_probes=designPrePost[(start<=designPrePost["start"]) & (end>=designPrePost["end"]) & (allchr==designPrePost["chr"])]

    if dataset['array_type'][index] == 'qChip 400K':
        num_probes=design400[(start<=design400["start"]) & (end>=design400["end"]) & (allchr==design400["chr"])]

    if dataset['array_type'][index] == 'qChip 1M':
        num_probes=design1M[(start<=design1M["start"]) & (end>=design1M["end"]) & (allchr==design1M["chr"])]

    if dataset['array_type'][index] == 'qChip 180K': 
        num_probes=design180[(start<=design180["start"]) & (end>=design180["end"]) & (allchr==design180["chr"])]

    if dataset['array_type'][index] == 'qChip Exon 8x60K (085824)': 
        num_probes=design60_24[(start<=design60_24["start"]) & (end>=design60_24["end"]) & (allchr==design60_24["chr"])]

    if dataset['array_type'][index] == 'qChip Exon 8x60K (085337)': 
        num_probes=design60_37[(start<=design60_37["start"]) & (end>=design60_37["end"]) & (allchr==design60_37["chr"])]


    if len(num_probes)>0:
        #adding the number of probes
        dataset['PROBES'][index] = str(len(num_probes))
        #assigning the coverage to True if the CNV has a minimum of three probes
        if len(num_probes) >= 3:
            dataset['COVERAGE'][index] = 'True'

    

#filtering only those CNVs covered
datasetCOVERED = pd.DataFrame(dataset, columns=dataset.columns)
datasetCOVERED = datasetCOVERED[~datasetCOVERED.COVERAGE.str.contains('False')]
datasetCOVERED.reset_index(inplace=True)


dataset.to_csv('{}/datasetCNV_VALIDATION_v2.csv'.format(os.getcwd()), header=True, sep="\t")
datasetCOVERED.to_csv('{}/datasetCNV_VALIDATION_COVERED.csv'.format(os.getcwd()), header=True, sep="\t")


# generating individual datasets with CNVs covered per sample

individual = pd.DataFrame(columns=datasetCOVERED.columns)
for sample in datasetCOVERED['Sample'].unique():
    individual=individual.append(datasetCOVERED[datasetCOVERED['Sample']==sample],ignore_index=True, sort=False)
    individual.to_csv('/Users/alba/Desktop/QGEN/cnv/EDcobertesIND/{}_EDcob.csv'.format(sample), header=True, sep="\t")
    individual = pd.DataFrame(columns=datasetCOVERED.columns)





#    C O V E R E D   C N V S   F R O M    A R R A Y S   -------------------------------------------------------------


color_print("### GENERATING THE WHOLE DATASET ###", color='red')
y=1
# creating an empty dataset where all CNVs from arrays files will be added
column_names = ['chrom','band','start','end','size','probes','type','genes','CNV %','Comments','Sample']
datasetCNVsARRAYS = pd.DataFrame(columns=column_names)

#Adding arrays CNVs data
for path in alldataarray:
    color_print('Adding ARRAY data {}     ({}/{})'.format(path.rsplit('/')[-1].rsplit('_')[2], y, len(alldataarray)), color='cyan')
    variable_length=len(datasetCNVsARRAYS)
    datasetCNVsARRAYS = datasetCNVsARRAYS.append(pd.read_csv(str(path), sep="\t"),ignore_index=True, sort=False)    #raises and error if the file is open!! BE CAREFUL!
    datasetCNVsARRAYS.loc[variable_length:len(datasetCNVsARRAYS),['Sample']] = path.rsplit('/')[-1].rsplit('_')[2]
    y+=1

datasetCNVsARRAYS.to_csv('{}/datasetCNVsARRAYS.csv'.format(os.getcwd()), header=True, sep="\t")


#creating three new columns to keep track of the number of exons and the coverage of each CNV
datasetCNVsARRAYS['numexons'] = [0] * len(datasetCNVsARRAYS) 
datasetCNVsARRAYS['COVERAGE'] = 'False'


#analysing each CNV individually
for i, cnva in datasetCNVsARRAYS.iterrows():
    chro = cnva['chrom']
    start = cnva['start']
    end = cnva['end']
    typee = cnva['type']

    #extracting the number of exons the CNV contains from the exoma dataset
    numexons = capexoma[(start<=capexoma["start"]) & (end>=capexoma["end"]) & (chro==capexoma["chr"])]

    if len(numexons)>0:
        #adding the number of exons and assigning the coverage to True if it has at least one exon
        datasetCNVsARRAYS['numexons'][i] = str(len(numexons))
        datasetCNVsARRAYS['COVERAGE'][i] = 'True'


#filtering only those CNVs covered
datasetCNVsARRAYSCOVERED = pd.DataFrame(datasetCNVsARRAYS, columns=datasetCNVsARRAYS.columns)
datasetCNVsARRAYSCOVERED = datasetCNVsARRAYSCOVERED[~datasetCNVsARRAYSCOVERED.COVERAGE.str.contains('False')]
datasetCNVsARRAYSCOVERED.reset_index(inplace=True)


datasetCNVsARRAYS.to_csv('{}/datasetCNVsARRAYS_v2.csv'.format(os.getcwd()), header=True, sep="\t")
datasetCNVsARRAYSCOVERED.to_csv('{}/datasetCNVsARRAYS_COVERED.csv'.format(os.getcwd()), header=True, sep="\t")





#    V A L I D A T I O N   E X O M E   D E P T H   A N D   A R R A Y S   --------------------------------------------


# using only the CNVs variants coming from ExomeDepth and arrays sufficiently covered
datasetCNVsARRAYS_COVERED = pd.read_csv('{}/datasetCNVsARRAYS_COVERED.csv'.format(os.getcwd()), sep="\t")
datasetCOB = pd.read_csv('{}/datasetCNV_VALIDATION_COVERED.csv'.format(os.getcwd()), sep="\t")

datasetCOB.drop(['Unnamed: 0'], axis=1, inplace=True)
datasetCNVsARRAYS_COVERED.drop(['Unnamed: 0'], axis=1, inplace=True)

truePdata = pd.DataFrame(columns = datasetCOB.columns)
falseNdata = pd.DataFrame(columns = ['chrom','band','start','end','size','probes','type','genes','CNV %','Comments'])
falsePdata = pd.DataFrame(columns = datasetCOB.columns)
segduppseudodata = pd.DataFrame(columns = datasetCOB.columns)

# matching names:  duplication == gain  and  deletion == loss
datasetCNVsARRAYS_COVERED['type'][datasetCNVsARRAYS_COVERED.type.str.contains('gain')] = 'duplication'
datasetCNVsARRAYS_COVERED['type'][datasetCNVsARRAYS_COVERED.type.str.contains('loss')] = 'deletion'

a=0
#analysing each CNV coming from ARRAYS individually
for ia, arrayy in datasetCNVsARRAYS_COVERED.iterrows():
    a+=1

    startarrayy = arrayy['start']
    endarrayy = arrayy['end']
    chrarrayy = arrayy['chrom']
    typearrayy = arrayy['type']  
    samplearrayy = arrayy['Sample']

    ARRAYinED1 = datasetCOB[(samplearrayy==datasetCOB['Sample']) & (startarrayy<=datasetCOB['start']) & (endarrayy>=datasetCOB['end']) & (typearrayy==datasetCOB['type']) & (chrarrayy[-1]==datasetCOB['chr'])]
    ARRAYinED2 = datasetCOB[((samplearrayy==datasetCOB['Sample']) & (startarrayy<datasetCOB['start']) & (endarrayy>=datasetCOB['start']) & (endarrayy<=datasetCOB['end']) & (typearrayy==datasetCOB['type']) & (chrarrayy[-1]==datasetCOB['chr']))]
    ARRAYinED3 = datasetCOB[((samplearrayy==datasetCOB['Sample']) & (startarrayy<=datasetCOB['end']) & (startarrayy>=datasetCOB['start']) & (endarrayy>datasetCOB['end']) & (typearrayy==datasetCOB['type']) & (chrarrayy[-1]==datasetCOB['chr']))]


    if (len(ARRAYinED1) == 0) and (len(ARRAYinED2) == 0) and (len(ARRAYinED3) == 0):
        #false negative data corresponds to those Exome Depth CNVs that are not predicted, but they are for arrays
        falseNdata = falseNdata.append(arrayy)





y=0
#analysing each CNV coming from EXOMEDEPTH individually
for i, cnvcob in datasetCOB.iterrows():
    y+=1
    idcnv = cnvcob['Sample']
    startcnv = cnvcob['start']
    endcnv = cnvcob['end']
    chrcnv = 'chr'+str(cnvcob['chr'])
    typecnv = cnvcob['type']   
    samplecnv = cnvcob['Sample']

    EDinARRAY1 = datasetCNVsARRAYS_COVERED[(samplecnv==datasetCNVsARRAYS_COVERED['Sample']) & (startcnv>=datasetCNVsARRAYS_COVERED['start']) & (endcnv<=datasetCNVsARRAYS_COVERED['end']) & (chrcnv==datasetCNVsARRAYS_COVERED['chrom']) & (typecnv==datasetCNVsARRAYS_COVERED['type'])]
    EDinARRAY2 = datasetCNVsARRAYS_COVERED[((samplecnv==datasetCNVsARRAYS_COVERED['Sample']) & (startcnv<datasetCNVsARRAYS_COVERED['start']) & (endcnv>=datasetCNVsARRAYS_COVERED['start']) & (endcnv<=datasetCNVsARRAYS_COVERED['end']) & (chrcnv==datasetCNVsARRAYS_COVERED['chrom']) & (typecnv==datasetCNVsARRAYS_COVERED['type']))]
    EDinARRAY3 = datasetCNVsARRAYS_COVERED[((samplecnv==datasetCNVsARRAYS_COVERED['Sample']) & (startcnv<=datasetCNVsARRAYS_COVERED['end']) & (startcnv>=datasetCNVsARRAYS_COVERED['start']) & (endcnv>datasetCNVsARRAYS_COVERED['end']) & (chrcnv==datasetCNVsARRAYS_COVERED['chrom']) & (typecnv==datasetCNVsARRAYS_COVERED['type']))]

    #SEGMENTAL DUPLICATIONS:  filtering CNVs falling within segmental duplications    
    sd1 = segdups[(startcnv>=segdups["chromStart"]) & (endcnv<=segdups["chromEnd"]) & (chrcnv==segdups["chrom"])]
    '''
    sd2 = segdups[(start<segdups["chromStart"]) & (end>segdups["chromEnd"]) & (allchr==segdups["chrom"])]
    sd3 = segdups[((start<segdups["chromStart"]) & (end>=segdups["chromStart"]) & (end<=segdups["chromEnd"]) & (allchr==segdups["chrom"]))]
    sd4 = segdups[((start<=segdups["chromEnd"]) & (start>=segdups["chromStart"]) & (end>segdups["chromEnd"]) & (allchr==segdups["chrom"]))]
    '''
    #considering pseudogenes
    pseudo=pseudogenes[(pseudogenes['cdsStart']>=startcnv) & (pseudogenes['cdsEnd']<=endcnv) & (pseudogenes['#chrom']==chrcnv)]
    pseudo.reset_index(inplace=True)

    #CNVs containing a superdup or a pseudogene
    if (len(sd1) > 0) or (len(pseudo) > 0): 
        segduppseudodata = segduppseudodata.append(cnvcob)

        
    else:   
        if (len(EDinARRAY1) > 0) or (len(EDinARRAY2) > 0) or (len(EDinARRAY3) > 0):
            #true positive data corresponds to those CNVs predicted by ExomeDepth and by arrays (coincide)
            truePdata = truePdata.append(cnvcob)

        else:
            #false positive data corresponds to those CNVs predicted by ExomeDepth but not by arrays (not coincide)
            falsePdata = falsePdata.append(cnvcob)






color_print('TRUE POSITIVE:     {}'.format(len(truePdata)), color='green')
color_print('FALSE POSITIVE:     {}'.format(len(falsePdata)), color='red')
color_print('FALSE NEGATIVE:     {}'.format(len(falseNdata)), color='magenta')
color_print('SEGDUPS AND PSEUDOGENES:     {}'.format(len(segduppseudodata)), color='white')

print('Total number of CNVs from arrays:', ia)
print('Total number of CNVs from Exome Depth with 3 at least:', y)


truePdata.to_csv('{}/truePdata.csv'.format(os.getcwd()), header=True, sep="\t")
falsePdata.to_csv('{}/falsePdata.csv'.format(os.getcwd()), header=True, sep="\t")
falseNdata.to_csv('{}/falseNdata.csv'.format(os.getcwd()), header=True, sep="\t")
segduppseudodata.to_csv('{}/segduppseudodata.csv'.format(os.getcwd()), header=True, sep="\t")

print()
print()
print()






stoptime = timeit.default_timer()
print('\n\nTime:',stoptime - starttime)




