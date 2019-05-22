import csv
import os
import datetime
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

#1. read blast result to identify shared gene cluster

#count pair of gene cluster in blast result
cluster_dict = {}                               #dictionary to keep pair of gene cluster as key and number of gene match as value
blast40 = []                                    #list to keep blast result with 40% identity

with open('./program/log/blast.tsv') as tsvfile:    #open and read blast result
  reader = csv.reader(tsvfile, delimiter='\t')

  query = ''
  subject = ''
  
  for row in reader:
      
    if float(row[2]) < 40:                          #filter out row[2] (percent identity) less than 40
        continue  
    elif (query, subject) == (row[0], row[1]):      #filter out pair of gene that match more than one time
        continue
    else:
        blast40.append(row)                         
        
        query = row[0]                              #gene of query sequence (main)
        query_cluster = query.split(':')[1]         #gene cluster noumber of guery
        subject = row[1]                            #gene of subject sequence (related)
        subject_cluster = subject.split(':')[1]     #gene cluster number of subject
        
        cluster_pair = query_cluster + '|' + subject_cluster
        #ex. Cordyceps_militaris_CM01_(ascomycetes)1|Cordyceps_confragosa_RCEF_1005_(ascomycetes)8
        
        try:
            cluster_dict[cluster_pair] += 1         #if cluster_pair already a key in cluster_dict, add value by 1
        except KeyError:
            cluster_dict[cluster_pair] = 1          #if cluster_pair not a key in cluster_dict, add cluster_pair as a key with value of 1

#write blast result with 40 %identity
with open('./program/log/blast40.tsv', 'wt') as blast_40_out_file:
    tsv_writer = csv.writer(blast_40_out_file, delimiter='\t')
    
    for r in blast40:
        tsv_writer.writerow(r)

#find number of genes in each gene cluster of main genome
main_cluster = {}
main_genecluster_fasta = list(SeqIO.parse("./program/log/main_genecluster.fasta", "fasta"))

for gene in main_genecluster_fasta:
    gene_cluster = gene.id.split(":")[1]
    try:
        main_cluster[gene_cluster] += 1         #if gene_cluster is already a key in main_cluster, add value by 1
    except KeyError:
        main_cluster[gene_cluster] = 1          #if gene_cluster not a key in main_cluster, add gene_cluster as a key with value of 1
            
#filter only pair that have more than 0.8 gene content similarity
shared_gene_cluster = []                            #list for shared gene cluster (>0.8 gene content)
all_group = []                                      #all group of shared gene cluster
group = []                                          #group of shared gene cluster (same main gene cluster)

for pair in list(cluster_dict.keys()):
    
    shared_gene = cluster_dict[pair]                #number of shared gene in pair of gene cluster
    main_genecluster, other_genecluster = pair.split("|")   
    n_main_gene = main_cluster[main_genecluster]    #number of gene in main gene cluster
    
    if (shared_gene / n_main_gene) > 0.8:           #more than 80% shared gene content
        shared_gene_cluster.append(pair)
        
        if main_genecluster not in group:           #new main gene cluster, add previous group to all_group
            all_group.append(group)
            group = []
            group.append(main_genecluster)          #add main gene cluster to group
            group.append(other_genecluster)         #add other gene cluster to group
            
        else:                                       #same main gene cluster as previous
            group.append(other_genecluster)         #add other gene cluster to group

all_group.append(group)                             #add last group to all_group
all_group.remove([])                                #remove empty list (from first group)

#create output folder with datetime as name
currentdatetime = datetime.datetime.now().strftime("%Y-%m-%d(%H_%M_%S)")
output_path = './output/' + currentdatetime
os.makedirs(output_path)

#write pair of gene cluster to file
with open(output_path + '/shared_gene_cluster_file.txt', 'w') as shared_file:
    for pair in shared_gene_cluster:
        shared_file.write(pair)
        shared_file.write('\n')

#write all group of gene cluster to file
with open('./program/log/group_gene_cluster_file.txt', 'w') as group_file:
    group_file.write(output_path)
    group_file.write('\n')
    
    for group in all_group:
        for cluster in group:
            group_file.write(cluster)
            group_file.write('\t')
        group_file.write('\n')

#create stacked bar graph

#find list of all organism
org_list = []
for group in all_group:
	for i in group:
		x = list(i.split('_('))[0]
		if x not in org_list:
			org_list.append(x)

import matplotlib.pyplot as plt
import pandas as pd

plt.close('all')

df = pd.DataFrame(data = None, index = org_list, columns = range(len(all_group)))

for c,group in enumerate(all_group):
	for i in group:
		x = list(i.split('_('))[0]
		df[c][x] = 1
		
df['SUM'] = df.count(axis='columns')
df.sort_values("SUM", ascending = True, inplace = True)
df.drop("SUM", axis=1, inplace = True)

ax = df.plot.barh(stacked=True,  legend = False, figsize = (10,5),
                  colormap = "nipy_spectral" )

ax.set_ylabel("Species")
ax.set_xlabel("Number of shared gene cluster")
plt.tight_layout()
plt.savefig(output_path + '/shared_gene_cluster_barplot.png')

