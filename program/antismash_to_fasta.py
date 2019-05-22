import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import xlrd
from Bio import Entrez

Entrez.email = "jennyjenny8@hotmail.com"

def find_org_name(assembly_acc):
    handle = Entrez.esearch(db="assembly", term = assembly_acc, retmode="xml")
    record = Entrez.read(handle)
    Id = record["IdList"][0]
    handle1 = Entrez.esummary(db="assembly", id=Id, report="full")
    record1 = Entrez.read(handle1)
    organism = record1['DocumentSummarySet']['DocumentSummary'][0]['Organism'].replace(' ', '_')
    return organism

    
#extract sequence of gene clusters from protein sequence by using antismash result
def antismash_genecluster(genome_sequence, antismash_result, org_name):
    n = antismash_result.nrows
    allgene = []
    
    for i in range(1,n):
        name = org_name + str(i)
        genes = list(antismash_result.cell(i,4).value.split(";"))
        
        for gene in genes:
            for CDS in genome_sequence:
                if gene in CDS.name:
                    a = SeqRecord(CDS.seq, id = CDS.id + ":" + name, description = '', name = '')
                    allgene.append(a)
    return allgene

#read input protein sequence (.fasta) and antismash result (.xls) of the main organism
for filename in os.listdir('./input/main'):                     #every file in the directory
    filename_w_ext = os.path.basename(filename)
    filename, file_extension = os.path.splitext(filename_w_ext) #separate filename and file extension

    if file_extension == ".faa":
        main_genome_fasta = list(SeqIO.parse( "./input/main/" + filename_w_ext, "fasta"))
        
        main_acc = "_".join(filename.split("_")[:2])
        organism = find_org_name(main_acc)

    if file_extension == ".xls":
        main_antismash = xlrd.open_workbook("./input/main/" + filename_w_ext).sheet_by_index(0)


main_genecluster_fasta = antismash_genecluster(main_genome_fasta, main_antismash, organism)
        
        

#read input protein sequence (.fasta) and antismash result (.xls) of the related organism
related_genecluster_fasta = []

for folder in os.listdir("./input/related"): 
    for filename in os.listdir("./input/related/" + folder):        #every file in the directory
        filename_w_ext = os.path.basename(filename)
        filename, file_extension = os.path.splitext(filename_w_ext) #separate filename and file extension

        if file_extension == ".faa":
            genome_fasta = list(SeqIO.parse( "./input/related/" + folder + '/' + filename_w_ext, "fasta"))

            acc = "_".join(filename.split("_")[:2])
            organism = find_org_name(acc)

        if file_extension == ".xls":
            antismash = xlrd.open_workbook("./input/related/" + folder + '/' + filename_w_ext).sheet_by_index(0)
            
    genecluster_fasta = antismash_genecluster(genome_fasta, antismash, organism)
    related_genecluster_fasta = related_genecluster_fasta + genecluster_fasta

#write fasta file of main and related gene cluster
SeqIO.write(main_genecluster_fasta,"./program/log/main_genecluster.fasta","fasta")
SeqIO.write(related_genecluster_fasta,"./program/log/related_genecluster.fasta","fasta")


