#coding:utf-8


import sys
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
from sklearn.preprocessing import MinMaxScaler




##### Input and ouput

inpfile=sys.argv[1]  ### input file, gene-SNP pairs
modulefile=sys.argv[2]  ### module file, conducted by Export_moduleinfo.py
bipfile=sys.argv[3]  ### module file, conducted by Export_moduleinfo.py
mod=sys.argv[4]  ### P/OD
oupfile=sys.argv[5]  ### onput file name


print('input: '+inpfile+'\t'+modulefile+'\t'+bipfile+'\n')
print('output: '+oupfile+'\n')

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
mre.columns=['gene','module','moduleset']  
ml = list(set(mre.module))
L2=mre

nl=GP.nodes()
Ng=set(mre.gene.astype('str'))  
N=len(Ng)

if mod=='P':
	print('mode: P-value')
	genetotal=[]
	pvaltotal=[]	
	genesettotal=[]


	genesets=[ list(L2[L2.module==i].gene.astype('str')) for i in ml]   

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


elif mod=='OD':
	print('mode: odd ratio')
	genetotal=[]
	odtotal=[]	
	genesettotal=[]


	genesets=[ list(L2[L2.module==i].gene.astype('str')) for i in ml]   

	num=1
	for gn in genesets:  
		g=set(gn).intersection(Ng)     
		M=len(g)                  
		for i in nl:	    
			nei=list(GP.neighbors(i))   
			neiF=set(nei).intersection(Ng)    
			n=len(neiF)     
			k=len(set(neiF).intersection(g))  
			od, pval = fisher_exact([[M-k,k],[N-M-n+k,n-k]])
			genetotal.append(i)
			odtotal.append(od)
			genesettotal.append(num)
		num=num+1
		print('1 finish')


	dfkda=pd.DataFrame({'gene':genetotal,'geneset':genesettotal,'od':odtotal})
	dfkda = dfkda.replace([np.inf, -np.inf], np.nan)
	dfkda = dfkda.replace(np.nan, 0)

	###### Normalized score
	gs=list(set(dfkda.geneset))
	ResTotal=pd.DataFrame()
	for i in gs:
		modres=dfkda[dfkda.geneset==i]
		modres['score']=modres['od']
		from sklearn.preprocessing import MinMaxScaler
		scaler = MinMaxScaler() 
		rownum='norm'+str(i)
		modres['norm'+str(i)]=scaler.fit_transform(modres['score'])
		if ResTotal.empty:
			ResTotal=modres.loc[:, ['gene', rownum]]
		else:
			ResTotal=pd.merge(ResTotal,modres.loc[:, ['gene', rownum]],on='gene')


######

		
ResTotal.to_csv(oupfile+"-KDA-EScore.txt",index=False,sep='\t' )



#####   SScore
ngenesets=len(genesets)
resori=ResTotal
#print(resori)
dg=pd.read_table(bipfile,sep='\t',header=0) 
dg.rename(columns={"node":"gene"},inplace=True)
dg['gene'] = dg.gene.astype('str')
#dg['DG']=dg['norm']+1
resweight=pd.merge(resori,dg,how = 'left',on='gene')
resweight=	resweight.fillna(1)


#print(len(resweight.columns))
#print(resweight[list(range(1,ngenesets+1))])

resweight_col = [column for column in resweight[list(range(1,ngenesets+1))]]
resweight_filt = resweight[resweight_col].multiply(resweight['norm'], axis="index")

		
scaler = MinMaxScaler() 
resweight_score = pd.DataFrame(columns=['gene'])
resweight_score['gene'] = resweight['gene']
resweight_score['SScore']=resweight_filt.mean(1)
resweight_score=resweight_score.sort(["SScore"],ascending=False)
resweight_score['norm']=scaler.fit_transform(resweight_score['SScore'])
resweight_score.to_csv(oupfile+"-SScore.txt",index=False,sep='\t' )


