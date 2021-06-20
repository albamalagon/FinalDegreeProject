'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'
'''

# LIBRARIES  -------------------------------------------------------------------------------------------

from pandas_ods_reader import read_ods
import timeit
import argparse
import sys
import warnings
import random
import scipy
import os
import re 
import pydot
import glob
import xlrd
import numpy as np
import pandas as pd
import matplotlib  
matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt 
from lazyme.string import color_print
from sklearn.metrics import roc_curve, auc
from pandas import read_csv
from pandas.plotting import scatter_matrix
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None  # default='warn'
from sklearn.ensemble import RandomForestClassifier
import graphviz
from sklearn import tree
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import classification_report, confusion_matrix, multilabel_confusion_matrix
from sklearn.impute import SimpleImputer
from sklearn import metrics
from numpy import where
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.datasets import make_classification
from sklearn.metrics import plot_confusion_matrix, ConfusionMatrixDisplay
from sklearn.tree import export_graphviz
from sklearn.model_selection import StratifiedKFold
from mlxtend.feature_selection import SequentialFeatureSelector
from keras.models import Sequential
from keras.layers import Dense
from keras.optimizers import SGD
import scipy.sparse
from sklearn.model_selection import learning_curve, ShuffleSplit
from keras.wrappers.scikit_learn import KerasClassifier






# LOADING DATASETS  ------------------------------------------------------------------------------------

#read files containing 20 runs
OK = glob.glob(f"{os.getcwd()}/Runs_qCarrier_Revisats/*.qCarrier.ALLsamples.AGRUPANT_TAULES.qAnotat.OK.debayan*.ods")
GF = glob.glob(f"{os.getcwd()}/allGFfiles/*.genome.FINAL.txt")
#read files containing 40 runs
OK40 = glob.glob(f"{os.getcwd()}/Runs_qCarrier_RevisatsMORE/*.qCarrier.ALLsamples.AGRUPANT_TAULES.qAnotat.OK.debayan*.ods")
GF40 = glob.glob(f"{os.getcwd()}/allGFfilesMORE/*.genome.FINAL.txt")


#read all filenames EXTERNAL DATA with many runs
OKe = glob.glob(f"{os.getcwd()}/external_data/allOKed/*.qCarrier.ALLsamples.AGRUPANT_TAULES.qAnotat.OK.debayan*.ods")
GFe = glob.glob(f"{os.getcwd()}/external_data/allGFed/*.genome.FINAL.txt")
#read all filenames EXTERNAL DATA with less runs
OK3runs = glob.glob(f"{os.getcwd()}/external_dataLESS/allOKed/*.qCarrier.ALLsamples.AGRUPANT_TAULES.qAnotat.OK.debayan*.ods")
GF3runs = glob.glob(f"{os.getcwd()}/external_dataLESS/allGFed/*.genome.FINAL.txt")
#read all filenames EXTERNAL DATA with only one run
OKonerun = glob.glob(f"{os.getcwd()}/onerun27/OK/*.qCarrier.ALLsamples.AGRUPANT_TAULES.qAnotat.OK.debayan*.ods")
GFonerun = glob.glob(f"{os.getcwd()}/onerun27/GF/*.genome.FINAL.txt")





# DEFINING VARIABLES  ----------------------------------------------------------------------------------


test_size=0.2
cumulative_importance=0.99
num_feature_importance=15
n_iterations=3
threshold=0.99









