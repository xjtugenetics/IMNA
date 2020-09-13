#coding:utf-8


import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import operator
import math
import pandas as pd 
from  scipy.stats import fisher_exact
import math
import csv
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')

##### Input and ouput

inpfile1=sys.argv[1]  ### input file, gene-SNP pairs
inpfile2=sys.argv[2]  ### module file, conducted by Export_moduleinfo.py
oupfile=sys.argv[3]  ### onput file name


############# PPI+GIANT
r1=pd.read_table(inpfile1,sep='\t',header=0)  
r1=r1.loc[:,['gene','norm']]
r1['gene'] = r1.gene.astype('str')

r2=pd.read_table(inpfile2,sep='\t',header=0)
r2=r2.loc[:,['gene','norm']]
r2['gene'] = r2.gene.astype('str')

rt=pd.merge(r1,r2,on='gene',how='outer')
rt=rt.fillna(0)

###
rt['mean']=rt.mean(1)
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
rt=rt.sort(['mean'],ascending=False)
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
rt['compscore']=scaler.fit_transform(rt['mean'])
rt=rt.sort(['compscore'],ascending=False)
rt.to_csv(oupfile+"_Composite_score.txt",index=False,sep='\t' )

######################
