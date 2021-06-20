'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'

Script to validate CNV calls between ExomeDepth, arrays and the developed method (CNV_predictions.py)

USAGE: python3 validation.py
'''

from CNV_libraries_variables import *

starttime = timeit.default_timer()


#    R E A D I N G    F I L E S   ----------------------------------------------------------------------------------

#these are the CNV calls between Exome Depth and Arrays (exomedepth_and_arrays.py)
truePdata = pd.read_csv('{}/truePdata.csv'.format(os.getcwd()), sep="\t")
falsePdata = pd.read_csv('{}/falsePdata.csv'.format(os.getcwd()), sep="\t")
segduppseudodata = pd.read_csv('{}/segduppseudodata.csv'.format(os.getcwd()), sep="\t")

#these are the CNV calls after implementing our strategy (CNV_validation_method.py)
truePdataM = pd.read_csv('{}/truePdataMETHOD.csv'.format(os.getcwd()), sep="\t")
falsePdataM = pd.read_csv('{}/falsePdataMETHOD.csv'.format(os.getcwd()), sep="\t")
noinfoM = pd.read_csv('{}/NOINFOdataMETHOD.csv'.format(os.getcwd()), sep="\t")
segdup_pseudoM = pd.read_csv('{}/superdups_pseudoMETHOD.csv'.format(os.getcwd()), sep="\t")

segduppseudodata.drop(['Unnamed: 0'], axis=1, inplace=True)
truePdata.drop(['Unnamed: 0'], axis=1, inplace=True)
falsePdata.drop(['Unnamed: 0'], axis=1, inplace=True)
truePdataM.drop(['Unnamed: 0'], axis=1, inplace=True)
falsePdataM.drop(['Unnamed: 0'], axis=1, inplace=True)
noinfoM.drop(['Unnamed: 0'], axis=1, inplace=True)
segdup_pseudoM.drop(['Unnamed: 0'], axis=1, inplace=True)


#only considering deletions
truePdata_indexes = truePdata[(truePdata['type'] != 'deletion')]
falsePdata_indexes = falsePdata[(falsePdata['type'] != 'deletion')]
segduppseudodata_indexes = segduppseudodata[(segduppseudodata['type'] != 'deletion')]
truePdataM_indexes = truePdataM[(truePdataM['type'] != 'deletion')]
falsePdataM_indexes = falsePdataM[(falsePdataM['type'] != 'deletion')]
noinfoM_indexes = noinfoM[(noinfoM['type'] != 'deletion')]
segdup_pseudoM_indexes = segdup_pseudoM[(segdup_pseudoM['type'] != 'deletion')]
# Step 2 : remove deletion indexes
truePdata = truePdata.drop(truePdata_indexes.index, axis=0) 
falsePdata = falsePdata.drop(falsePdata_indexes.index, axis=0) 
segduppseudodata = segduppseudodata.drop(segduppseudodata_indexes.index, axis=0) 
truePdataM = truePdataM.drop(truePdataM_indexes.index, axis=0) 
falsePdataM = falsePdataM.drop(falsePdataM_indexes.index, axis=0) 
noinfoM = noinfoM.drop(noinfoM_indexes.index, axis=0) 
segdup_pseudoM = segdup_pseudoM.drop(segdup_pseudoM_indexes.index, axis=0) 




#printing heads of all datasets
color_print('Coincide ED and ARRAYS:', color='magenta', end='')
print(len(truePdata))
print(truePdata.head())

color_print('NO coincide ED and ARRAYS:', color='red', end='')
print(len(falsePdata))
print(falsePdata.head())

color_print('Compatible ED and METHOD:', color='green',end='')
print(len(truePdataM))
print(truePdataM.head())

color_print('NO compatible ED and METHOD:', color='cyan',end='')
print(len(falsePdataM))
print(falsePdataM.head())

color_print('NO info ED and METHOD:', color='yellow',end='')
print(len(noinfoM))
print(noinfoM.head())

color_print('SegDups and pseudogenes arrays:', color='blue',end='')
print(len(segduppseudodata))
print(segduppseudodata.head())

color_print('SegDups and pseudogenes METHOD:', color='white',end='')
print(len(segdup_pseudoM))
print(segdup_pseudoM.head())


#merging the coincident calls between ExomeDepth and arrays with each of the possibilities:

#with compatibility cases between ExomeDepth and the method
COINCIDEIXENtp_OK = truePdata.merge(truePdataM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENtp_OK.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENtp_OK.csv', header=True, sep="\t")
#with INcompatibility cases between ExomeDepth and the method
COINCIDEIXENtp_KO = truePdata.merge(falsePdataM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENtp_KO.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENtp_KO.csv', header=True, sep="\t")
#with no information cases between ExomeDepth and the method
COINCIDEIXENtp_NI = truePdata.merge(noinfoM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENtp_NI.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENtp_NI.csv', header=True, sep="\t")
#with segmental duplication (or pseudogenes) cases between ExomeDepth and the method
COINCIDEIXENtp_SD = truePdata.merge(segdup_pseudoM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENtp_SD.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENtp_SD.csv', header=True, sep="\t")


#merging the NO coincident calls between ExomeDepth and arrays with each of the possibilities:

#with compatibility cases between ExomeDepth and the method
COINCIDEIXENfp_OK = falsePdata.merge(truePdataM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner',indicator=True)
COINCIDEIXENfp_OK.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENfp_OK.csv', header=True, sep="\t")
#with INcompatibility cases between ExomeDepth and the method
COINCIDEIXENfp_KO = falsePdata.merge(falsePdataM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner',indicator=True)
COINCIDEIXENfp_KO.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENfp_KO.csv', header=True, sep="\t")
#with no information cases between ExomeDepth and the method
COINCIDEIXENfp_NI = falsePdata.merge(noinfoM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner',indicator=True)
COINCIDEIXENfp_NI.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENfp_NI.csv', header=True, sep="\t")
#with segmental duplication (or pseudogenes) cases between ExomeDepth and the method
COINCIDEIXENfp_SD = falsePdata.merge(segdup_pseudoM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner',indicator=True)
COINCIDEIXENfp_SD.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENfp_SD.csv', header=True, sep="\t")


#merging the segmental duplication calls between ExomeDepth and arrays with each of the possibilities:

#with segmental duplication (or pseudogenes) cases between ExomeDepth and the method
COINCIDEIXENsd_SD = segduppseudodata.merge(segdup_pseudoM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENsd_SD.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENsd_SD.csv', header=True, sep="\t")
#with no information cases between ExomeDepth and the method
COINCIDEIXENsd_ni = segduppseudodata.merge(noinfoM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENsd_ni.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENsd_ni.csv', header=True, sep="\t")
#with compatibility cases between ExomeDepth and the method
COINCIDEIXENsd_tp = segduppseudodata.merge(truePdataM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENsd_tp.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENsd_tp.csv', header=True, sep="\t")
#with INcompatibility cases between ExomeDepth and the method
COINCIDEIXENsd_fp = segduppseudodata.merge(falsePdataM, on=['chr', 'start', 'end', 'width', 'type', 'nexons', 'BF',
       'reads.expected', 'reads.observed', 'reads.ratio', 'SegDup',
       'exons.hg19', 'RefSeqGenes', 'Sample', 'PROBES', 'COVERAGE',
       'array_type'], how = 'inner' ,indicator=True)
COINCIDEIXENsd_fp.to_csv('/Users/alba/Desktop/QGEN/cnv/COINCIDEIXENsd_fp.csv', header=True, sep="\t")



color_print('\n \n \n------------------------------------------------------------------------------------------------------------------------------------------------------------\n \n \n'.format(), color='red')



color_print('Coincide ED and ARRAY + compatible METHOD:                 --> TP: {}'.format(len(COINCIDEIXENtp_OK)), color='magenta')
print(COINCIDEIXENtp_OK.head())
color_print('Coincide ED and ARRAY + NO compatible METHOD:              --> FP: {}'.format(len(COINCIDEIXENtp_KO)), color='magenta')
print(COINCIDEIXENtp_KO.head())
color_print('Coincide ED and ARRAY + NO INFO METHOD:     {}'.format(len(COINCIDEIXENtp_NI)), color='magenta')
print(COINCIDEIXENtp_NI.head())




color_print('No coincide ED and ARRAY + compatible METHOD:                 --> FP: {}'.format(len(COINCIDEIXENfp_OK)), color='cyan')
print(COINCIDEIXENfp_OK.head())
color_print('No coincide ED and ARRAY + NO compatible METHOD:              --> TP: {}'.format(len(COINCIDEIXENfp_KO)), color='cyan')
print(COINCIDEIXENfp_KO.head())
color_print('No coincide ED and ARRAY + NO INFO METHOD:     {}'.format(len(COINCIDEIXENfp_NI)), color='cyan')
print(COINCIDEIXENfp_NI.head())

color_print('SEGDUP or PSEUDOGENE METHOD:     {}'.format(len(COINCIDEIXENsd_SD)), color='cyan')
print(COINCIDEIXENsd_SD.head())


print('\n')
print('Total number of CNVs from arrays: 247')
print('Total number of CNVs from Exome Depth covered: 992') #992 o 570

color_print('\n \n \n------------------------------------------------------------------------------------------------------------------------------------------------------------\n \n \n'.format(), color='red')




t = PrettyTable(["call coincident exome depth and d'arrays", 'method OK', 'method KO', 'method NI', 'method SD'])
t.add_row([len(truePdata), len(COINCIDEIXENtp_OK), len(COINCIDEIXENtp_KO), len(COINCIDEIXENtp_NI), len(COINCIDEIXENtp_SD)])
print(t)

print('\n\n')

t2 = PrettyTable(["call NO coincident exome depth and d'arrays", 'method OK', 'method KO', 'method NI', 'method SD'])
t2.add_row([len(falsePdata), len(COINCIDEIXENfp_OK), len(COINCIDEIXENfp_KO), len(COINCIDEIXENfp_NI), len(COINCIDEIXENfp_SD)])
print(t2)

print('\n\n')
t3 = PrettyTable(["total segmental duplications", 'method SD', 'method NI', 'method OK', 'method KO'])
t3.add_row([len(segduppseudodata), len(COINCIDEIXENsd_SD), len(COINCIDEIXENsd_ni), len(COINCIDEIXENsd_tp), len(COINCIDEIXENsd_fp) ])
print(t3)

print('\n')


color_print('\nONLY CONSIDERING DELETIONS:\n'.format(), color='magenta')




color_print('Total TP:     {}'.format(len(COINCIDEIXENtp_OK)+len(COINCIDEIXENfp_KO)), color='green')
color_print('Total FP:     {}'.format(len(COINCIDEIXENtp_KO)+len(COINCIDEIXENfp_OK)), color='red')
color_print('Total NI:     {}'.format(len(COINCIDEIXENtp_NI)+len(COINCIDEIXENfp_NI)), color='yellow')
color_print('Total SD:     {}'.format(len(COINCIDEIXENtp_SD)+len(COINCIDEIXENfp_SD)+len(COINCIDEIXENsd_SD)), color='white')





stoptime = timeit.default_timer()
print('\n\nTime:',stoptime - starttime)
