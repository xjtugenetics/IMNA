#coding:utf-8

###################### import
import os
import re
import sys
import pandas as pd

inpmoduledir=sys.argv[1]
oup=sys.argv[2]


def file_name(file_dir):    ######
	for root, dirs, files in os.walk(file_dir):
		return(files)
		
flist=file_name(inpmoduledir)

fn=[]

for f in flist:
	fname=re.split('\W|_',f)[0]
	fn.append(fname)

ftype=set(fn)	
moduleinfo=pd.DataFrame()
modname=[]
n1=1
n3=1
n2=101
for fm in ftype:
	fileset=[file.find(fm) for file in flist]
	fileindex=[i for i,v in enumerate(fileset) if v==0]
	for file in fileindex:
		f1=open(inpmoduledir+'/'+flist[file],'r').readlines()
		genelist=[gene.strip() for gene in f1]
		modulelist=[n2]*len(genelist)
		modulesetlist=[n3]*len(genelist)
		df1=pd.DataFrame({'moduleset':modulesetlist,'module':modulelist,'gene':genelist})
		moduleinfo=pd.concat([moduleinfo,df1])
		modname.append(flist[file])
		n2=n2+1
	n3=n3+1
modulen=list(set(moduleinfo.module))
df2=pd.DataFrame({'module':modulen,'name':modname})	
	
moduleinfo.to_csv(oup+'.txt', sep='\t', header=True, index=False)	
df2.to_csv(oup+'info.txt', sep='\t', header=True, index=False)

