'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'
'''

# LIBRARIES  -------------------------------------------------------------------------------------------

import timeit
import sys
import warnings
import os
import re 
import glob
import numpy as np
import pandas as pd
import matplotlib  
matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt 
from lazyme.string import color_print
from scipy.stats import binom 
pd.options.mode.chained_assignment = None
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
warnings.simplefilter("ignore", category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
from prettytable import PrettyTable


# LOADING DATASETS  ------------------------------------------------------------------------------------

GF = glob.glob(f"{os.getcwd()}/run596/Resultats_SNVs_Indels/*.genome.FINAL.txt") 
ED = glob.glob(f"{os.getcwd()}/run596/Resultats_CNVs_ExomeDepth/*.ExomeDepth*.txt") 
RM = glob.glob(f"{os.getcwd()}/repeatmasker_17032021_modified.txt") #using USCS GRCH37 
SD = glob.glob(f"{os.getcwd()}/segmentaldups.txt") #using USCS GRCH37 
pseudogenes = glob.glob(f"{os.getcwd()}/pseudo.txt") #using USCS GRCH37 


# DEFINING VARIABLES  ----------------------------------------------------------------------------------

#'min_variants_affected' is the number that defines the minimum number of variants a CNV must contain to have significant results
min_variants_affected = 3
#'min_variants_control' is the minimum number of control variables to make a control distribution
min_variants_control = 3
#'min_prob_freq_random' defines the minimum probability to have a variant by chance
min_prob_freq_random = 0.000001


# VALIDATION DATASETS  --------------------------------------------------------------------------------

GFval = glob.glob(f"{os.getcwd()}/GF_validation/*.genome.FINAL.txt") 
EDval = glob.glob(f"{os.getcwd()}/ED_validation/*_Solapaments_*.txt")
EDvalcob =glob.glob(f"{os.getcwd()}/EDcobertesIND/*_EDcob.csv") 
design = glob.glob(f"{os.getcwd()}/design/*.bed") 
array = glob.glob(f"{os.getcwd()}/arrays_alba/*_Selection.xls") 
exoma = glob.glob(f"{os.getcwd()}/20210127-jrodriguez-mostres_amb_array_i_exoma.csv") 
captura_ex = glob.glob(f"{os.getcwd()}/MedExome_hg19_capture_targets.bed") 


# data obtained from the CNV_validation_method.py
TPdata = glob.glob(f"{os.getcwd()}/truePdata.csv") 
FPdata = glob.glob(f"{os.getcwd()}/truePdata.csv") 
