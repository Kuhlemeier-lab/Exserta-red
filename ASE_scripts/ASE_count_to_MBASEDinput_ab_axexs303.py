###### THIS IS AN INTERACTIVE SCRIPT I RUN ON LINUX WITH PYTHON3 #####
#### THIS SCRIPT WAS MODIFIED FROM AN ORIGINAL SCRIPT WRITTEN BY M. MOSER, PLEASE SEE https://doi.org/10.1016/j.cub.2018.10.019 ####
# JUPYTER NOTEBOOK OPTIONAL #

#BEFORE RUNNING THIS SCRIPT
#you will need your counts files and your genome annotation file (gff) and python3
#then followed the instructions here to make sure that the jupyter notebook could run ipython with a python3 kernel
#https://ipython.readthedocs.io/en/latest/install/kernel_install.html
	#python3 -m pip install ipykernel
	#python3 -m ipykernel install --user
	
# parse ASE read counts for analysis with MBASED: 
#Take only first 8 columns of file as rest is empty or unimportant to this analysis
#from your terminal:
cut -f 1-8 F1s.strict.counts > F1s.strict.counts.short
cut -f 1-8 4.strict.counts > axex1.strict.counts.short
cut -f 1-8 5.strict.counts > axex2.strict.counts.short
cut -f 1-8 6.strict.counts > axex3.strict.counts.short


#pip3 install pandas
#pip3 install matplotlib
#pip3 install scipy


python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas import *
import collections
from collections import defaultdict
from __future__ import division
from collections import Counter
from scipy.stats.mstats import chisquare


## Synchronize the ASE.counts.short files:
ASErep1 = "axex1.strict.counts.short"
ASErep2 = "axex2.strict.counts.short"
ASErep3 = "axex3.strict.counts.short"


def loadASE(filename): 
    file_dict = {}
    with open(filename, "r") as f: 
        next(f)  #skip the first line
        count = 0
        for l in f:
            count += 1
            l = l.strip()
            l = l.split('\t')
            head = '_'.join(l[:5])
            rest = l[5:8] #is now set to include 3 replicates !!! has to be changed for more or less replicates         
            file_dict[head] = rest
            #print file_dict
           
    return file_dict 
    
#%%time
#sync the file for each position

t1 = loadASE(ASErep1)
t2 = loadASE(ASErep2)
t3 = loadASE(ASErep3)

#python3 specific
allpos = set(list(t1.keys()) + list(t2.keys()) + list(t3.keys()))

print(len(t1), len(t2), len(t3))


print(len(allpos))

out = []



e = 0
for i in allpos: 
    e += 1
    head = i.split('_')[0]
    pos = i.split('_')[1]
    if i in t1 or i in t2 or i in t3: 
        #if i not in defaultdict, defaultdict gets extended by item i with default value 0
        out.append([i, t1.get(i, ['0', '0', '0']), t2.get(i, ['0', '0', '0']), t3.get(i, ['0', '0', '0'])])
        
print(len(out))
print(out)

print(out[1])

print('\t'.join(out[1][0].split('_')[0:2] + out[1][0].split('_')[3:5]))

combined_counts = []
with open('F1s.strict.counts.short', "w") as f: 
		e = 0
		for i in out: 
			e += 1
			#print i
			id_loc = '\t'.join(i[0].split('_')[0:2] + i[0].split('_')[3:5])
			counts = '\t'.join(['\t'.join(i[1]), '\t'.join(i[2]) ,'\t'.join(i[3])])
			f.write(id_loc +'\t'+ counts + '\n')
			#print id_loc, counts
			li = id_loc.split('\t')
			li.extend(counts.split('\t'))
			combined_counts.append(li)
			if e < 10: 
				print(id_loc, counts)
				print(li)

#IMPORTANT list to use in the next step when looking up genes matching the posistions
combined_counts[0]

#load genome annotation of interest

#check for correct version of annotation, I prefer to sort my gff before loading it
gene_list = []
with open("peaxi162AQ_PeaxHIC303.cds.swapphase.sorted.AB.fixedtabs.gff") as gff:
	ei = 0
	for i in gff: 
		if i.startswith("##"): 
			continue
		ei += 1
		#if ei < 10: 
		#    print i
		i = i.strip()
		i = i.split("\t")
		if i[2] != "gene":
			continue
		#gene = []
		#gene.append(i[0])
		#gene.append(i[3:5])
		gene_name = i[8].split(";")[0][3:]
		gene_list.append([gene_name, i[0], int(i[3]), int(i[4])])
     
        
gene_df = pd.DataFrame(gene_list, columns=['gene', 'sequence', 'first', 'last'])
print(gene_df.shape)
print(gene_df.head(10))
print(gene_df.loc[gene_df["sequence"]=="Peax302Scf85298"]) #Enter a scaffold here to test


%%time

#write to output

#what do positions look like?
c = 0
k = 0
with open('F1s.strict.counts.short', 'w') as outf: 
    for e, i in enumerate(combined_counts):
        #print e, i
        #only SNPs within gene annotations get printed to file
        result = gene_df.loc[(gene_df["sequence"] == i[0])&(gene_df["first"] <= int(i[1]))&(gene_df["last"] >= int(i[1]))]
        rest = 0
        #print(gene_df.loc[gene_df['gene'] == result])
        if result.shape[0] > 0:
            c += 1
            stringo = '\t'.join(i)
            res = '\t'.join([str(result.values[0][i]) for i in (0,2,3)])
            fin = res + '\t' + stringo + '\n'
            #print stringo, res
            outf.write(fin)
            if c % 10000 == 0:
                k += 1
                print(str(k) + ' * 10K processed...')
        
print(c)


#afterwards, sort with "$bedtools sort -i yourcountsfile.counts" to not mix up SNPs when treated in R
