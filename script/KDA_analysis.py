#coding:utf-8

##### Input and ouput

inpfile=sys.argv[1]  ### input file, gene-SNP pairs
modulefile=sys.argv[2]  ### module file, conducted by Export_moduleinfo.py
oupfile=sys.argv[3]  ### onput file name


import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import operator
import math
import pandas as pd 
import random
from  scipy.stats import fisher_exact
import math
import csv
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')

######## Define function
def NetConstructWDG(file1):
	DGP = nx.DiGraph()
	#f1.readline()
	for line in f1:
		l=line.strip().split()
		if len(l)==3 and l[0] !=l[1]:
			DGP.add_edge(l[0], l[1], weight=float(l[2]))		
	return(DGP)

def NetWDGtoWG(DGP):	
	GP=DGP.to_undirected()
	return(GP)
	
############### Open interaction file	
f1=open(inpfile)
DGP=NetConstructWDG(f1)
GP=NetWDGtoWG(DGP)


############### Open module file
mre=pd.read_table(modulefile,sep='\t',header=0) 
ml = list(set(mre.module))
L2=mre

nl=GP.nodes()
Ng=set(mre.gene).intersection(set(GP.nodes()))  
N=len(Ng)

genetotal=[]
pvaltotal=[]	
genesettotal=[]


genesets=[ list(L2[L2.module==i].gene) for i in ml]
num=1
for gn in genesets:
	g=set(gn).intersection(Ng)
	M=len(g)
	for i in nl:	
		nei=list(GP.neighbors(i))
		neiF=set(nei).intersection(Ng)
		n=len(neiF)
		k=len(set(neiF).intersection(g))
		pval = fisher_exact([[M-k,k],[N-M-n+k,n-k]])[1]
		genetotal.append(i)
		pvaltotal.append(pval)
		genesettotal.append(num)
	num=num+1
	print('1 finish')


dfkda=pd.DataFrame({'gene':genetotal,'geneset':genesettotal,'pval':pvaltotal})
dfkda.to_csv(oupfile+"_KDA-module.txt",index=False,sep='\t' )	 ##模块基因集

###### Normalized score
gs=list(set(dfkda.geneset))
ResTotal=pd.DataFrame()
for i in gs:
	modres=dfkda[dfkda.geneset==i]
	modres['score']=-(modres['pval'].apply(np.log10))
	from sklearn.preprocessing import MinMaxScaler
	scaler = MinMaxScaler() 
	rownum='norm'+str(i)
	modres['norm'+str(i)]=scaler.fit_transform(modres['score'])
	if ResTotal.empty:
		ResTotal=modres.loc[:, ['gene', rownum]]
	else:
		ResTotal=pd.merge(ResTotal,modres.loc[:, ['gene', rownum]],on='gene')

		
testscore=	ResTotal.mean(1)
ResTotal['score']=ResTotal.mean(1)
ResTotal=ResTotal.sort(["score"],ascending=False)
ResTotal['norm']=scaler.fit_transform(ResTotal['score'])
ResTotal.to_csv(oupfile+"_KDA-score.txt",index=False,sep='\t' )


