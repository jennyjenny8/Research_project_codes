python ./program/antismash_to_fasta.py

makeblastdb.exe -in ./program/log/related_genecluster.fasta -dbtype prot

blastp.exe -query ./program/log/main_genecluster.fasta -db ./program/log/related_genecluster.fasta -evalue 0.00001 -outfmt 6 -out ./program/log/blast.tsv

python ./program/shared_gene.py

python ./program/gene_cluster_diagram.py

python ./program/phylogenetic_tree.py

PAUSE