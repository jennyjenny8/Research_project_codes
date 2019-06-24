#3. Draw gene cluster diagram for each group of shared gene cluster
#==================================================================
import csv
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
#==================================================================

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

#open genbank file of gene cluster
def open_genbank_cluster(gene_cluster, path):       #ex. gene_cluster = Cordyceps_militaris_CM01_(ascomycetes)1
    cluster_no = gene_cluster.split(")")[1]         #cluster_no = 1
    length_no = len(cluster_no)                     
    if length_no == 1:
        cluster_no = "00" + cluster_no              #add 00 in front of cluster_no
    elif length_no == 2:
        cluster_no = "0" + cluster_no               #add 0 in front of cluster_no

    for file in os.listdir(path):                   #read one file in the specific folder
        genome_id = file.split("cluster")[0]        #ex. file = JH126399.1.cluster001.gbk, genome_id = JH126399.1.
        break

    cluster_file = path + '/' + genome_id + 'cluster' + cluster_no + ".gbk" #filename of gene cluster number ...
    record = SeqIO.read(cluster_file, "genbank")    #read genbank file of gene cluster
    return record                                   #output record from function


for group in all_group:                             #each group of shared gene cluster
    main = group[0]                                 #main gene cluster in group

    q = ''
    group_blast = []
    
    for n in range(len(group)-1):                   #add empty list to group_blast for each related gene cluster
        group_blast.append([])

    main_blast = []
    
    for row in blast40:
        #ex. row = [EGX95479.1:Cordyceps_militaris_CM01_(ascomycetes)1,
        #           OAA74936.1:Cordyceps_confragosa_RCEF_1005_(ascomycetes)76,
        #           51.703,323,121,6,4,298,16,331,8.75e-099,294]
        
        query_cluster = row[0].split(':')[1]        #Cordyceps_militaris_CM01_(ascomycetes)1
        subject_cluster = row[1].split(':')[1]      #Cordyceps_confragosa_RCEF_1005_(ascomycetes)76
        
        if query_cluster == main:
            q = row[0]                              #found main 
            if subject_cluster in group:            #subject gene cluster is in group of shared gene cluster
                pos_in_group = group.index(subject_cluster) #find position in group
                blast = [row[0].split(':')[0], row[1].split(':')[0], row[2]]    #blast = [query gene, subject gene, percent identity]
                group_blast[pos_in_group -1].append(blast)  #blast to group_blast at position found -1 (-1 because no main in group)

                query_gene = row[0].split(':')[0]   
                if query_gene not in main_blast:    #add query gene to main_blast
                    main_blast.append(query_gene)

        else:
            if q != '' :                            #already found main, now query_cluster not main so break from for loop (for row in blast40:)
                break
    
    n_main_blast = len(main_blast)                  #number of gene in main cluster that match with other
            
    records = dict()        
    records["A_record"] = open_genbank_cluster(main, "./input/main/cluster") #open genbank file of main gene cluster
    
    #open genbank file of other gene cluster in group
    records_list = ["A_record"]
    alphabet_list = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    i = 1
    for cluster in group[1:]:                       #related gene cluster in group
        foldername = cluster[0] + "_" + cluster.split("_")[1]   #ex. C_confragosa
        name = alphabet_list[i] + "_record"
        records_list.append(name)
        #open genbank file of related gene cluster and keep in records
        records[name] = open_genbank_cluster(cluster, "./input/related/" + foldername + "/cluster")
        i += 1
    
    #genome diagram
    gd_diagram = GenomeDiagram.Diagram(main)
    feature_sets = {}
    max_len = 0

    main_record_color = dict()                      #color of gene in main gene cluster

    n_group = len(group)                            #number of gene cluster in group

    l = 0
    for record_name in records_list:                #each gene cluster
        record = records[record_name]               #record = opened genbank file of gene cluster
        max_len = max(max_len, len(record))

        if n_group == 3:                              
            max_track = 4                           #to expand track
        elif n_group == 2:
            max_track = 3
        else:
            max_track = n_group -1
            
        gd_track_for_features = gd_diagram.new_track(max_track -l,  #create new track
                                name = record.name, height = 0.4,
                                start = 0, end = len(record))
        
        gd_feature_set = gd_track_for_features.new_set()

        c = 0
        m = 0
        for feature in record.features:
            if feature.type != "CDS":               #Exclude feature that is no CDS
                continue

            name = feature.qualifiers["protein_id"][0]  #get protein id of feature

            if l == 0:                              #main gene cluster
                if name in main_blast:              #feature in main that match with other
                    color_value = m/n_main_blast    #rank of feature in main blast hit

                    #color range [red,orange,yellow,green,cyan,blue,purple,magenta,pink]
                    if color_value < 0.125:
                        color = colors.linearlyInterpolatedColor(colors.red, colors.orange,0, 0.125, color_value)
                    elif color_value < 0.25:
                        color = colors.linearlyInterpolatedColor(colors.orange, colors.yellow,0.125, 0.25, color_value)    
                    elif color_value < 0.375:
                        color = colors.linearlyInterpolatedColor(colors.yellow, colors.green,0.25, 0.375, color_value)
                    elif color_value < 0.5:
                        color = colors.linearlyInterpolatedColor(colors.green, colors.cyan,0.375, 0.5, color_value)
                    elif color_value < 0.625:
                        color = colors.linearlyInterpolatedColor(colors.cyan, colors.blue,0.5, 0.625, color_value)
                    elif color_value < 0.75:
                        color = colors.linearlyInterpolatedColor(colors.blue, colors.purple,0.625, 0.75, color_value)
                    elif color_value < 0.875:
                        color = colors.linearlyInterpolatedColor(colors.purple, colors.magenta,0.75, 0.875, color_value)
                    else:
                        color = colors.linearlyInterpolatedColor(colors.magenta, colors.pink,0.875, 1, color_value)
                    main_record_color[name] = color
                    m += 1
                    
                else:                               #feature in main that don't match with other
                    color = colors.black
                    main_record_color[name] = color
            else:                                   #other gene cluster
                inpair = False
                for pair in group_blast[l-1]:       #find if this feature is in blast hit
                    if name in pair:
                        p = pair
                        inpair = True
                        break
                    
                if inpair == False:                 #feature not in blast hit
                    color = colors.black
                else:
                    try:                            #set color same with main hit
                        color = main_record_color[p[0]]
                    except KeyError:
                        color = colors.white        #if blast hit not in feature of main
                    
            #set start, end and strand of feature
            f = SeqFeature(FeatureLocation(feature.location.nofuzzy_start, feature.location.nofuzzy_end), strand=feature.strand)
            #add feature to diagram
            gd_feature_set.add_feature(f, sigil = "BIGARROW", color = color,
                                       arrowshaft_height = 1.0 , label=True, name =  name)

            c += 1
        l += 1

    gd_diagram.draw(format = "linear", pagesize = 'A4', fragments = 1 ,start = 0, end = max_len)    #draq diagram

    gd_diagram.write(output_path + '/' + main + '.pdf', "PDF")  #write PDF file
