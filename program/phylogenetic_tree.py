#4. Phylogenetic analysis of shared gene cluster
#================================================

#change the location of mafft.bat file below
mafft_exe = "C:/Users/User/Desktop/tools/mafft-win/mafft.bat"

#===========================================================
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import csv

from Bio.Align.Applications import MafftCommandline
from Bio.Alphabet import IUPAC
from Bio.Nexus import Nexus

#open gene cluster group file

group_file = open('./program/log/group_gene_cluster_file.txt').read()

lines = list(group_file.split('\n'))

output_path = lines[0]

all_group = []

for line in lines[1:]:
    group = []
    list_group = list(line.split('\t'))
    list_group.remove('')
    
    for cluster in list_group:
        group.append(cluster)

    all_group.append(group)

all_group.remove([])

#open blast40
blast40_file = open('./program/log/blast40.tsv')    #open and read blast result
blast40 = list(csv.reader(blast40_file, delimiter='\t'))

while blast40.count([]) > 0:
    blast40.remove([])

#gene cluster fasta file
main_fasta = list(SeqIO.parse("./program/log/main_genecluster.fasta", "fasta"))
related_fasta = list(SeqIO.parse("./program/log/related_genecluster.fasta", "fasta"))


nexus_filename = []

for group in all_group:
    group_fasta = []
    main = group[0]
    main_short = main[0]  + list(main.split("_"))[1] + list(main.split(")"))[1]
    main_name = main[0] + list(main.split("_"))[1]

    main_cluster_seq_list = []
    for sequence in main_fasta:
        if main == list(sequence.name.split(":"))[1]:
             main_cluster_seq_list.append(str(sequence.seq))
    main_cluster_seq = ''.join(main_cluster_seq_list)

    group_main_fasta = SeqRecord(Seq(main_cluster_seq),id = main_name, name = '', description = '')
    group_fasta.append(group_main_fasta)

    for cluster in group[1:]:
        cluster_name = cluster[0] + '' + list(cluster.split("_"))[1]

        #search for this pair in blast result
        pair_blast = []
        for row in blast40:
            if (list(row[0].split(":"))[1] == main) and (list(row[1].split(":"))[1] == cluster):
                pair_blast.append([list(row[0].split(":"))[0], list(row[1].split(":"))[0]])

        #check order of blast hit
        reverse = False
        if (float(pair_blast[0][1][3:]) - float(pair_blast[1][1][3:])) > 0:
            reverse = True
        
        cluster_seq_list = []
        for sequence in related_fasta:
            if cluster == list(sequence.name.split(":"))[1]:
                cluster_seq_list.append(str(sequence.seq))

        if reverse == False:    
            cluster_seq = ''.join(cluster_seq_list)
        else:
            cluster_seq_list.reverse()                  #reverse order of genes
            cluster_seq = ''.join(cluster_seq_list)

        fasta = SeqRecord(Seq(cluster_seq),id = cluster_name, name = '', description = '')
        group_fasta.append(fasta)

    SeqIO.write(group_fasta, "./program/log/" + main_short + ".fasta", 'fasta')

    #align group fasta using mafft
    mafft_cline = MafftCommandline(mafft_exe, input = "./program/log/" + main_short + ".fasta")
    stdout, stderr = mafft_cline()
    aligned_file_name = "./program/log/" + main_short + "_aligned.fasta"
    
    with open(aligned_file_name, "w") as handle:
        handle.write(stdout)
    
    AlignIO.convert(aligned_file_name,'fasta', aligned_file_name.replace(".fasta",".nex"),'nexus',alphabet=IUPAC.protein)
    nexus_filename.append(aligned_file_name.replace(".fasta",".nex"))

#combind nexus files
nexi =  [(fname, Nexus.Nexus(fname)) for fname in nexus_filename]

combined = Nexus.combine(nexi)
nexi_filename = open('./program/log/COMBINED.nex', 'w')
combined.write_nexus_data(nexi_filename)
nexi_filename.close()

#convert nexus into phylip
AlignIO.convert('./program/log/COMBINED.nex', "nexus", "./program/log/COMBINED.phy", "phylip")

#Build neighbor joining tree
aln = AlignIO.read('./program/log/COMBINED.phy', 'phylip')
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(aln)
Phylo.write(tree,output_path + '/COMBINED_nj.nhx', 'newick')

