# IMNA
###  1. Introduction
------------
</br>

> **`IMNA`** is a integrative multi-omics network-based approach to capture genetic-driven regulatory networks for human complex diseases.  
- This method can combine functional data from multiple biological scales to understand molecular mechanisms of disease and identify potential key genes.  
- This pipeline provide several scripts facilitating data access, integration and analysis.

</br>

###  2. Download and Configure
------------
</br><b>2.1 Download</b>
<pre>
	git clone https://github.com/xjtugenetics/IMNA.git
</pre>

</br><b>2.2 Configure</b>
<pre>
	export IMNA_tk=/path/to/IMNA
</pre>

</br>

###  3. Tutorial
------------
</br>  

- **Workflow**  

![IMNA workflow](https://github.com/xjtugenetics/IMNA/blob/master/workflow.png)

</br><b>3.1 Export gene set and module information </b>
```bash
Script:
		1-Export_module.py
		
Usage:
		python3 ${IMNA_tk}/script/1-Export_module.py <inpmoduledir> <oup>

```

##### Details for parameters:  
<pre>
inpmoduledir          Directory of all gene module (gene list)
oup                   Prefix of output file
</pre>

##### Output file(s):  
> ***oup.txt***: Gene set (module)  
> ***oupinfo.txt***: Module information  

</br><b>3.2 Constract bipartite based on SNP-gene pairs</b>
```bash
Script:
		2-Constract_bipartite.py

Usage:
		python3 ${IMNA_tk}/script/2-Constract_bipartite.py <GS pair> <oup>

```

##### Details for parameters:
<pre>
GS pair               Gene-SNP pairs 
oup                   Prefix of output file
</pre>

##### Output file(s):
> ***oup_DG_snp.txt***: Raw SNP degree score  
> ***oup_DG_gene.txt***: Raw gene degree score  
> ***oup_DG_gene.nor.txt***: Normalized gene degree score  

</br><b>3.3 Key driver analysis and Signature score analysis</b>
```bash
Script:
		3-Enricment_combine_SScore.py

Usage:
		python3 ${IMNA_tk}/script/3-Enricment_combine_SScore.py <inter_score> <modulefile> <gene_norm_dg> <mod> <oup>

```

##### Details for parameters:
<pre>
inter_score           Interaction score of database (PPI, GIANT)
modulefile            module file, conducted by 1-Export_module.py
gene_norm_dg          Normalized gene degree score, conducted by 2-Constract_bipartite.py
mod                   P-value or odd ratio	  
oup                   Prefix of output file
</pre>

##### Output file(s):
> ***oup-KDA-EScore.txt***: Gene enrichment score  
> ***oup-SScore.txt***: Combined gene signiture score 


</br><b>3.4 CS - Composite score analysis </b>
```bash
Script:
		4-Composite_score.py

Usage:
		python3 ${IMNA_tk}/script/5-Composite_score.py <PPI_SScore> <GIANT_SScore> <oup>

```

##### Details for parameters:
<pre>
PPI_SScore            Gene signiture score of PPI network (step 3 output)
GIANT_SScore          Gene signiture score of GIANT network (step 4 output)
oup                   Prefix of output file
</pre>

##### Output file(s):
> ***oup_Composite_score.txt***: Gene composite score (PPI and GIANT network)  

</br>  


</br><b>3.5 optional - Gene ID conversion (between symbol and entrz id)  </b>
```bash
Script:
		convert_geneid.r

Usage:
		Rscript ${IMNA_tk}/script/convert_geneid.r <inpfile> <oup> <column> <gene id type>

```

##### Details for parameters:
<pre>
inpfile           	  File before conversion
oup          		  Output file name
column           	  The number of target column
gene id           	  type SYMBOL or ENTREZID
</pre>

##### Output file(s):
> ***oup***: Gene id after conversion

</br>  


###  4. Example
------------
</br>

* Calculate gene composition score of 6 gene set modlue accoding to gene interaction from PPI and GIANT netwrok

```sh
###Step1:  

		python3 ${IMNA_tk}/script/1-Export_module.py ${IMNA_tk}/data/geneset module

		## "module.txt"
		gene	module	moduleset
		ADAMDEC1	101	1
		ASRGL1	101	1
		C10orf116	101	1
		CD55	101	1
		COL1A2	101	1
		COL3A1	101	1
		CTHRC1	101	1
		CXCL9	101	1
		CXCL10	101	1

		## "moduleinfo.txt"
		module	name
		101	39-signatures-ori.txt
		102	14-signatures-ori.txt
		103	300-signatures-ori.txt
		104	263-signatures-ori.txt
		105	126-signature-ori.txt
		106	100-signatures-ori.txt

```
```sh
###(optional)  

		Rscript ${IMNA_tk}/script/convert_geneid.r  module.txt  module.txt.corv 1 SYMBOL

		## "module.txt"
		gene	module	moduleset
		ADAMDEC1	101	1
		ASRGL1	101	1
		C10orf116	101	1
		CD55	101	1
		COL1A2	101	1
		COL3A1	101	1
		CTHRC1	101	1
		CXCL9	101	1
		CXCL10	101	1

		## "module.txt.corv"
		ENTREZID	module	moduleset
		13	106	6
		10157	104	4
		23460	104	4
		22	103	3
		8714	103	3
		3983	103	3
		52	103	3
		56	103	3
		59	103	3

```
```sh
###Step2:  

		python3 ${IMNA_tk}/script/2-Constract_bipartite.py snp.gene.pairs bip

		## "snp.gene.pairs"
		ENTREZID	SNP
		79575	rs13343778
		79575	rs113634768
		79575	rs11540855
		79575	rs34546260
		79575	rs62126253
		79575	rs55924783
		79575	rs4808616
		79575	rs28473003
		79575	rs12982058

		## "bip_DG_snp.txt"
		DG	node
		0.006993006993006993	chr17:43989159:I
		0.006993006993006993	rs112058513
		0.006993006993006993	chr17:46685076:I
		0.006993006993006993	rs1457919
		0.006993006993006993	rs62070651
		0.006993006993006993	rs7207826
		0.006993006993006993	rs138922609
		0.006993006993006993	rs16940671
		0.006993006993006993	rs4748775

		## "bip_DG_gene.txt"
		DG	node
		0.00022660321776569228	580
		0.00022660321776569228	752
		0.00022660321776569228	148345
		0.00022660321776569228	246744
		0.00022660321776569228	1021
		0.00022660321776569228	1267
		0.1907999093587129	1394
		0.00045320643553138455	1545
		0.00022660321776569228	83468

		## "bip_DG_gene.nor.csv"
		DG	node	norm
		0.22003172445048721	284058	2.0
		0.1907999093587129	1394	1.8670103092783505
		0.1470654883299343	4137	1.6680412371134021
		0.10151824155903014	9884	1.4608247422680412
		0.10061182868796738	474170	1.4567010309278352
		0.09811919329254476	51326	1.445360824742268
		0.05030591434398369	389170	1.2278350515463918
		0.04146838885112169	8631	1.1876288659793814
		0.024699750736460458	100506084	1.111340206185567

```
```sh
###Step3:  

		python3 ${IMNA_tk}/script/3-Enricment_combine_SScore.py PPI.filter.txt   module.txt bip_DG_gene.nor.txt P PPI
		python3 ${IMNA_tk}/script/3-Enricment_combine_SScore.py GIANT.filter.txt  module.txt  bip_DG_gene.nor.txt P GIANT

		## "PPI-KDA-EScore.txt"
		gene	norm1	norm2	norm3	norm4	norm5	norm6
		23136	0.0	0.0	0.0	0.0	0.0	0.0
		55586	0.0	0.0	0.0	0.0	0.0	0.0
		4958	0.0	0.0	0.0	0.0	0.0	0.0
		9793	0.0	0.0	0.0	0.0	0.0	0.0
		8543	0.0	0.0	0.0	0.0	0.0	0.0
		117	0.0	0.0	0.0	0.0	0.0	0.0
		200845	0.0	0.0	0.0	0.0	0.0	0.28086621546950413
		79139	0.0	0.0	0.0	0.0	0.0	0.0
		3912	0.0	0.49879745804727266	0.0	0.0	0.0	0.0

		## "PPI-SScore.txt"
		gene	SScore	norm
		4137	0.47177454297544336	1.0
		51226	0.3023830979296318	0.6409483140453625
		3292	0.28409462722290135	0.6021830373278301
		5158	0.2638202235224249	0.559208264732837
		26276	0.26273712445041764	0.5569124666908819
		23331	0.23212051246337354	0.49201576456289586
		6155	0.2275463373024836	0.48232008422363687
		55839	0.22354248460319304	0.4738332916255484
		1021	0.2191097252448783	0.464437364218449

		## "GIANT-KDA-EScore.txt"
		gene	norm1	norm2	norm3	norm4	norm5	norm6
		100462799	0.0	0.0	0.0	0.0	0.0	0.05798346369741207
		402360	0.06340642377742223	0.0	0.0	0.0	0.0	0.0
		404266	0.0	0.0	0.0	0.0	0.0	0.0
		441067	0.0	0.0	0.0	0.0	0.0	0.0
		400860	0.0	0.1461858411609449	0.0	0.0	0.0	0.0
		100128252	0.0	0.0	0.0	0.0	0.0	0.0
		143188	0.0	0.0	0.0	0.0	0.0	0.0
		340107	0.0	0.0	0.0	0.0	0.0	0.0
		28982	0.0	0.0	0.0	0.0	0.0	0.0

		## "GIANT-SScore.txt"
		gene	SScore	norm
		10197	0.6911514953381394	1.0
		3093	0.6615120899310956	0.9571159064156506
		56942	0.649453035818157	0.939668133829933
		8678	0.5406093711987038	0.7821865030245152
		3612	0.5280639792253968	0.764035067256921
		55839	0.5148197866889375	0.7448725643530102
		10951	0.5049760346912975	0.7306300255405549
		9055	0.44075870302240827	0.6377164861761186
		130507	0.4081344528874548	0.5905137377844765
		

```
```sh
###Step4:  

		python3 ${IMNA_tk}/script/4-Composite_score.py PPI-SScore.txt GIANT-SScore.txt result

		## "result_Composite_score.txt"
		gene	norm_x	norm_y	mean	compscore
		10197	0.3688901262723868	1.0	0.6844450631361934	1.0
		55839	0.4738332916255484	0.7448725643530102	0.6093529279892793	0.8902875640554194
		4137	1.0	0.1532682979117107	0.5766341489558553	0.8424841963409919
		3093	0.14321740008504674	0.9571159064156506	0.5501666532503486	0.8038141888692014
		56942	0.05267380011244613	0.939668133829933	0.49617096697118956	0.7249244587983239
		8678	0.20451842336073286	0.7821865030245152	0.49335246319262405	0.720806518688345
		9055	0.3255048079610211	0.6377164861761186	0.4816106470685699	0.7036512833649254
		6155	0.4823200842236369	0.4469615605789528	0.4646408224012949	0.6788577307757409
		3612	0.15652586829854412	0.7640350672569209	0.4602804677777325	0.672487088545402



```
```sh
###(optional) 

		convert_geneid.r result_Composite_score.txt result_Composite_score.txt.corv 1 ENTREZID 

		## "result_Composite_score.txt.corv"
				SYMBOL	norm_x	norm_y	mean	compscore
		NAT2	0	0	0	0
		ADA	0	0.0139823333689054	0.00699116668445271	0.0102143576760106
		CDH2	0	0.00698273440278748	0.00349136720139374	0.00510101889755177
		AKT3	0	0.00698273440278748	0.00349136720139374	0.00510101889755177
		GAGE12F	0	0	0	0
		RNA18S5	0	0	0	0
		RNA28S5	0	0.0267595577506105	0.0133797788753052	0.0195483605565037
		TRF-GAA1-4	0	0.00417259333462038	0.00208629666731019	0.00304815795989612
		ANO1-AS2	0	0.0551490756412158	0.0275745378206079	0.0402874376714162

```
```sh
###(optional) 

		sort  -k5 -n -r  result_Composite_score.txt.corv  | head 

```

</br>


###  5. License
------------
</br>  

> This software is distributed under the terms of GPL 2.0
</br>


###  6. Source
------------
</br>  

> [https://github.com/xjtugenetics/IMNA](https://github.com/xjtugenetics/IMNA)
</br>


###  7. Contact
------------
</br>  

#### Author
> **Yi-Xiao Chen**, **Yu Rong**, **Feng Jiang**, **Jia-Bin Chen**, **Yuan-Yuan Duan**, **Shan-Shan Dong**, **Dong-Li Zhu**, **Hao Chen**, **Tie-Lin Yang**, **Zhijun Dai**, **Yan Guo**
> Key Laboratory of Biomedical Information & Genetics Center, School of Life Science and Technology, Xi'an Jiaotong University, Xi'an, Shaanxi Province, 710049, P. R. China  
> [:email:](guoyan253@mail.xjtu.edu.cn) guoyan253@mail.xjtu.edu.cn  
</br>


###  8. Maintainer
------------
</br>  

> **Yi-Xiao Chen**  
> You can contact [:email:](chenyixiao@stu.xjtu.edu.cn) chenyixiao@stu.xjtu.edu.cn
  when you have any questions, suggestions, comments, etc.
> Please describe in details, and attach your command line and log messages if possible.  
</br>


###  9. Requiremnets
------------
</br>  
	
- **Python** ( >= 3.2 )
	- numpy( >=1.10.4)
	- mpmath( >=0.19)
	- pandas( >=0.18.0)
	- scipy( >= 0.17.0)
	- rpy2( >= 2.9.0)
	- networkx( >= 2.2)
	- sklearn( >= 0.17.1)
	- sys
	- os

- **R** \( >= 3.3.2 \)
	- org.Hs.eg.db
	- clusterProfiler
</br>


###################### Thank you! #############################

