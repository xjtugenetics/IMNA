#coding:utf-8


#####################	
import sys
import os
import pandas as pd 
from  scipy.stats import fisher_exact
import math
import csv
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
import networkx as nx
from networkx.algorithms import bipartite 
from sklearn.preprocessing import scale
from networkx.algorithms import bipartite 

#####################

inpfile=sys.argv[1]  ### input file, gene-SNP pairs
oupprefix=sys.argv[2]  ### onput file name




#####################

def NetBIPconstruct(dataframe,top_nodes,bottom_nodes):	
	b = nx.Graph()
	gene=set(top_nodes)
	SNP=set(bottom_nodes)
	b.add_nodes_from(SNP, bipartite=0)
	b.add_nodes_from(gene, bipartite=1)
	kgp=[]
	for index, row in dataframe.iterrows():
		kgp.append((str(row[top_nodes]), row[bottom_nodes]))
	b.add_edges_from(kgp)	
	return(b)

#####################	Extract gene and SNP

SNP2genetotal=pd.read_table(inpfile,sep='\t',header=0)
SNP2genetotal.columns = ['gene', 'SNP']
SNP2gene=SNP2genetotal
b = nx.Graph()
gene=set(SNP2gene.gene)
snp=set(SNP2gene.SNP)

b.add_nodes_from(snp, bipartite=0)
b.add_nodes_from(gene, bipartite=1)
kgp=[]
for index, row in SNP2gene.iterrows():
    kgp.append((str(row["SNP"]), row["gene"]))

	
##################### Add edges
count_all=[ kgp.count(i) for i in kgp if kgp.count(i)]
count_list=set(count_all)

for n in count_list:	
	kgp_wei1= [ i for i in kgp if kgp.count(i) == n]
	b.add_edges_from(kgp_wei1,weight=n)
	
bottom_nodes = set(n for n,d in b.nodes(data=True) if d['bipartite']==1)
top_nodes = set(b) - bottom_nodes

bip_dg=bipartite.degree_centrality(b,bottom_nodes)	
	
node=[]
dgvalue=[]
for  i in bip_dg:
	node.append(i)
	dgvalue.append(bip_dg[i])
	
dg=pd.DataFrame({'node':node,'DG':dgvalue})
dg_snp=dg[dg.node.isin(top_nodes)]
dg_gene=dg[dg.node.isin(bottom_nodes)]

dg_snp.sort(['DG'])
dg_gene.sort(['DG'])
dg_snp.to_csv(oupprefix+"_DG_snp.txt",index=False,sep='\t' )
dg_gene.to_csv(oupprefix+"_DG_gene.txt",index=False,sep='\t' )

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler() 
dg_gene['norm']=scaler.fit_transform(dg_gene['DG'])+1
dg_gene=dg_gene.sort(['norm'],ascending=False)
dg_gene.to_csv(oupprefix+"_DG_gene.nor.txt",index=False,sep='\t' )




