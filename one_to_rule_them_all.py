import pandas as pd
from os import walk
import argparse
from collections import defaultdict
import networkx as nx
import matplotlib
import matplotlib.pyplot
import csv
import numpy as np

def main(args, fantom_files: list):

    results_df = pd.read_csv(args.results_file, delimiter=",",header=None) 
    paths_df = pd.read_csv(args.path_file, delimiter=",",header=None) # go and genes
    cond1 = pd.read_csv(f"condition_C9_vs_C9GC_unique_genes.csv", delimiter=",",header=None)
    cond2 = pd.read_csv(f"condition_C9KO_vs_C9GC_unique_genes.csv", delimiter=",",header=None)

    results_df = results_df.iloc[1:] #remove first row
    paths_df = paths_df.iloc[2:,1:]
    
    results_df = results_df.iloc[(-results_df[5].astype(float).abs()).argsort()] # sort by biggest nes
     
    separator: str = ''
    
    intersec = []
    onion = []
    
    for i in range(len(cond1.iloc[0])):
        for j in range(len(cond2.iloc[0])):
            if cond1.iloc[0,i] == cond2.iloc[0,j]:
                intersec.append(cond1.iloc[0,i])
            if (cond1.iloc[0,i] not in onion ):
                onion.append(cond1.iloc[0,i])

    genes_genes = pd.read_excel(f"SG_enriched_collection.xlsx", header=None).iloc[1:]
    gene_list = []
    prot_list = []

    gene_list = list(set(genes_genes[1]).intersection(set(onion)))
    prot_list = list(set(onion).difference(set(gene_list)))
    print(len(gene_list) + len(prot_list) == len(onion))
    ################## all rezults
    rez_dict: defaultdict= defaultdict(int)

    nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])
    quantity = 0
    separator = ''
    not_found = []
    for gene in onion:
        entered = 0
        for file in fantom_files:
            
            if gene in file:
                entered = 1
                if len(file) > 3:
                    i=0
                    graph_file = ''
                    for junc in file:
                        
                        if 0 != i % 2 or len(file) == i:
                            graph_file += junc + ","
                        else: 
                            graph_file += junc

                else:
                    graph_file = f"{args.phantom_folder}/{separator.join(file)}"
                graph_results = pd.read_csv(graph_file, delimiter=",", skiprows=1)

                for i in range(len(graph_results)):
                    if (graph_results.iloc[i, 3] >= .95 and not graph_results.iloc[i].isnull().iloc[-3]):
                        rez_dict[' '.join(sorted([gene, graph_results.iloc[i, -3]]))] += 1
                    else:
                        break
        if entered == 0:
            not_found.append(gene)

    nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])

    for key, it in rez_dict.items():
        wan, two = key.split(' ')
        nodes.loc[quantity] = [wan, two, 'Undirected', it]
        quantity += 1
 
   ######################### MY deg
   
    all_names: defaultdict= defaultdict(list)
    for i in range(len(nodes)):
        
        if (nodes.iloc[i,0] in onion):
            if(nodes.iloc[i,1] not in all_names[nodes.iloc[i,0]]):
                all_names[nodes.iloc[i,0]].append(nodes.iloc[i,1])
                
        elif (nodes.iloc[i,1] in onion):
            if(nodes.iloc[i,0] not in all_names[nodes.iloc[i,1]]):
                all_names[nodes.iloc[i,1]].append(nodes.iloc[i,0])
    
    sec_degree: defaultdict= defaultdict(list)

    for key1 in all_names.keys():
        for key2 in all_names.keys():
            if key1 != key2:
                if (list(set(all_names[key1]).intersection(set(all_names[key2]))) != [] and len(sec_degree[' '.join(sorted([key1, key2]))]) == 0):
                    sec_degree[' '.join(sorted([key1, key2]))].append(set(all_names[key1]).intersection(set(all_names[key2])))

    for key, it in sec_degree.items():
        print(f"{key}: {it}")

    third_degree: defaultdict= defaultdict(list)

    for key1 in sec_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(sec_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(third_degree[' '.join(sorted(names + [key2]))]) == 0):
                    third_degree[' '.join(sorted(names + [key2]))].append(sec_degree[key1][0].intersection(set(all_names[key2])))
    print('\n\nthird_degree')
    for key, it in third_degree.items():
        print(f"{key}: {it}")

    forth_degree: defaultdict= defaultdict(list)

    for key1 in third_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(third_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(forth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    forth_degree[' '.join(sorted(names + [key2]))].append(third_degree[key1][0].intersection(set(all_names[key2])))
    print('\n\nforth_degree')
    for key, it in forth_degree.items():
        print(f"{key}: {it}")

    fith_degree: defaultdict= defaultdict(list)

    for key1 in forth_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(forth_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(fith_degree[' '.join(sorted(names + [key2]))]) == 0):
                    fith_degree[' '.join(sorted(names + [key2]))].append(forth_degree[key1][0].intersection(set(all_names[key2])))
    print('\n\nfith_degree')
    for key, it in fith_degree.items():
        print(f"{key}: {it}")
    
    sixth_degree: defaultdict= defaultdict(list)

    for key1 in fith_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(fith_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(sixth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    sixth_degree[' '.join(sorted(names + [key2]))].append(fith_degree[key1][0].intersection(set(all_names[key2])))
    print('\n\nsixth_degree')
    for key, it in sixth_degree.items():
        print(f"{key}: {it}")
    
    seventh_degree: defaultdict= defaultdict(list)

    for key1 in sixth_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(sixth_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(seventh_degree[' '.join(sorted(names + [key2]))]) == 0):
                    seventh_degree[' '.join(sorted(names + [key2]))].append(sixth_degree[key1][0].intersection(set(all_names[key2])))
    print('\n\nseventh_degree')
    for key, it in seventh_degree.items():
        print(f"{key}: {it}")
    
    eighth_degree: defaultdict= defaultdict(list)

    for key1 in seventh_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(seventh_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(eighth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    eighth_degree[' '.join(sorted(names + [key2]))].append(seventh_degree[key1][0].intersection(set(all_names[key2])))
    print('\n\n eighth_degree')
    for key, it in eighth_degree.items():
        print(f"{key}: {it}")

    fig = matplotlib.pyplot.figure(figsize=(40,40))
    G = nx.from_pandas_edgelist(nodes, source='source', target='target', edge_attr='weight')
    nx.draw_kamada_kawai(G, ax=fig.add_subplot(), font_size=8, with_labels=True, node_size=100, edge_color='orange')
    matplotlib.use("Agg")
    fig.savefig('all_cond_graph' + '.png')
    
    nodes_soft = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])

    all_eigth = []
    for a in eighth_degree.keys():
        all_eigth = a.split(' ')
    for a in eighth_degree.values():
        all_eigth.extend(list(a[0]))

    colour = ["red", "blue", "green", "yellow"]
    cond1_set = set(cond1.iloc[0])
    cond2_set = set(cond2.iloc[0])
    cond1_int = set(all_eigth).intersection(cond1_set)
    cond2_int = set(all_eigth).intersection(cond2_set)
    quantity = 0
    for key, it in eighth_degree.items():
        for gene in key.split(' '):
            for b in list(it[0]):
                if gene in cond1_int and gene in cond2_int:
                    
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]

                elif gene in cond1_int:
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]
                
                elif gene in cond2_int:
                    
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]
                
                else:                
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]
                quantity += 1
        
    quantity = 0
    nodesss = []
    colours = []
    for unique in all_eigth:
        if unique in cond1_int and unique in cond2_int:
            nodesss.append(unique)
            colours.append(colour[0])

        elif unique in cond1_int:
            nodesss.append(unique)
            colours.append(colour[1])
        
        elif unique in cond2_int:
            
            nodesss.append(unique)
            colours.append(colour[2])
        
        else: 
            nodesss.append(unique)
            colours.append(colour[3])
        quantity += 1

    fig = matplotlib.pyplot.figure(figsize=(12,12))
    G = nx.from_pandas_edgelist(nodes_soft, edge_attr=True)
    nx.draw_kamada_kawai(G, ax=fig.add_subplot(), font_size=8, with_labels=True, nodelist=nodesss, node_size=200, edge_color='orange', node_color=colours)
    matplotlib.use("Agg")
    fig.savefig('eigth' + '.png')
    
    nodes_soft = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])

    all_seventh = []
    for a in seventh_degree.keys():
        for b in a.split(' '):
            if b not in all_seventh:
                all_seventh.append(b)
    for a in seventh_degree.values():
        for b in list(a[0]):
            if b not in all_seventh:
                all_seventh.append(b)

    colour = ["red", "blue", "green", "yellow"]
    cond1_set = set(cond1.iloc[0])
    cond2_set = set(cond2.iloc[0])
    cond1_int = set(all_seventh).intersection(cond1_set)
    cond2_int = set(all_seventh).intersection(cond2_set)
    quantity = 0
    weird_ones = []
    for key, it in seventh_degree.items():
        print(key.split(' '))
        for gene in key.split(' '):
            for b in list(it[0]):
                if gene in cond1_int and gene in cond2_int:
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]

                elif gene in cond1_int:
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]
                
                elif gene in cond2_int:
                    
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]
                
                else:          
                    nodes_soft.loc[quantity] = [gene, b, 'Undirected', 1]
                quantity += 1
                
    newdf = pd.DataFrame(columns=['', 'baseMean',	'log2FoldChange',	'lfcSE'	,'stat'	,'pvalue'	,'padj'	,'gene_name', 'cond'])
    
    cond1_excell = pd.read_csv(f"condition_C9_vs_C9GC_diff.csv", header=None).iloc[1:]
    cond2_excell = pd.read_csv(f"condition_C9KO_vs_C9GC_diff.csv", header=None).iloc[1:]
    quantity = 0
    for i in range(len(cond1_excell.iloc[:])):
        if cond1_excell.iloc[i][7] in weird_ones:
            newdf.loc[quantity] = list(cond1_excell.iloc[i]) + [1]
            quantity += 1

    for i in range(len(cond2_excell.iloc[:])):
        if cond2_excell.iloc[i][7] in weird_ones:
            newdf.loc[quantity] = list(cond2_excell.iloc[i]) + [2]
            quantity += 1 

    quantity = 0
    nodesss = []
    colours = []
    for unique in all_seventh:
        if unique in cond1_int and unique in cond2_int:
            nodesss.append(unique)
            colours.append(colour[0])

        elif unique in cond1_int:
            nodesss.append(unique)
            colours.append(colour[1])
        
        elif unique in cond2_int:
            
            nodesss.append(unique)
            colours.append(colour[2])
        
        else: 
            weird_ones.append(unique)   
            nodesss.append(unique)
            colours.append(colour[3])
        quantity += 1

    c1 = "RPS19 EIF3K RPS6 EIF3B EIF3F EIF3D PABPC1 RPS3 RPS24 EIF4B EIF1".split(' ')# : Green (found in both)
    c2 = 'EIF3H SLBP EIF3M EIF3I HNRNPD CNBP EIF4A1 EIF3L EIF4H CSDE1 EIF2A EIF3E HABP4 EIF4E3'.split(' ')#  : Red (found in first condition TDP43)

    c3 = 'EIF3C RPS11 CLNS1A EIF3G RAN POP7 SNRPF FBL RNASEL NCOA3 MSI1 MSI2 FABP5 EIF4G3 CALR SAMD4A TARDBP RBM47 CNOT2'.split(' ') #: Blue (in second condition)
    c4 = c1 + c2 + c3
    
    nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])
    mark = pd.read_csv('/home/rubs/Downloads/network2.csv', delimiter=",",header=None).iloc[1:]
    quantity = 0
    unique = []
    source = []

    for i in range(len(mark.iloc[:])):
        if mark.iloc[i,1] not in unique and mark.iloc[i,1] not in c4:
            source.append(mark.iloc[i,1])
            unique.append(mark.iloc[i,1])
        if mark.iloc[i,11] not in unique and mark.iloc[i,11] not in c4:
            unique.append(mark.iloc[i,11])

        nodes.loc[quantity] = [mark.iloc[i,1], mark.iloc[i,11], 'Undirected', 1]
        quantity += 1 
    c_all = (len(c1) * ['green']) + (len(c2) * ['red']) + (len(c3) * ['blue']) + (len(unique) * ['yellow'])

    fig = matplotlib.pyplot.figure(figsize=(40,40))
    G = nx.from_pandas_edgelist(nodes, edge_attr=True)
    nx.draw_kamada_kawai(G, ax=fig.add_subplot(), font_size=6, with_labels=True, nodelist=c4 + unique, node_size=150, edge_color='orange', node_color=c_all)
    matplotlib.use("Agg")
    fig.savefig('mark_wazowsky' + '.png')

    ############## mark deg

    all_names: defaultdict= defaultdict(list)
    
    for i in range(len(nodes)):
        
        if (nodes.iloc[i,0] in unique):
            if(nodes.iloc[i,1] not in all_names[nodes.iloc[i,0]]):
                all_names[nodes.iloc[i,0]].append(nodes.iloc[i,1])
                
        elif (nodes.iloc[i,1] in unique):
            if(nodes.iloc[i,0] not in all_names[nodes.iloc[i,1]]):
                all_names[nodes.iloc[i,1]].append(nodes.iloc[i,0])
    
    sec_degree: defaultdict= defaultdict(set)

    for key1 in all_names.keys():
        for key2 in all_names.keys():
            if key1 != key2:
                if (list(set(all_names[key1]).intersection(set(all_names[key2]))) != [] and len(sec_degree[' '.join(sorted([key1, key2]))]) == 0):
                    sec_degree[' '.join(sorted([key1, key2]))] = set(all_names[key1]).intersection(set(all_names[key2]))


    third_degree: defaultdict= defaultdict(list)

    for key1 in sec_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (sec_degree[key1].intersection(set(all_names[key2])) != [] and len(third_degree[' '.join(sorted(names + [key2]))]) == 0):
                    third_degree[' '.join(sorted(names + [key2]))] = sec_degree[key1].intersection(set(all_names[key2]))


    forth_degree: defaultdict= defaultdict(list)

    for key1 in third_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(third_degree[key1].intersection(set(all_names[key2]))) != [] and len(forth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    forth_degree[' '.join(sorted(names + [key2]))] = third_degree[key1].intersection(set(all_names[key2]))


    fith_degree: defaultdict= defaultdict(list)

    for key1 in forth_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(forth_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(fith_degree[' '.join(sorted(names + [key2]))]) == 0):
                    fith_degree[' '.join(sorted(names + [key2]))].append(forth_degree[key1][0].intersection(set(all_names[key2])))

    
    sixth_degree: defaultdict= defaultdict(list)

    for key1 in fith_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(fith_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(sixth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    sixth_degree[' '.join(sorted(names + [key2]))].append(fith_degree[key1][0].intersection(set(all_names[key2])))

    
    seventh_degree: defaultdict= defaultdict(list)

    for key1 in sixth_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(sixth_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(seventh_degree[' '.join(sorted(names + [key2]))]) == 0):
                    seventh_degree[' '.join(sorted(names + [key2]))].append(sixth_degree[key1][0].intersection(set(all_names[key2])))

    
    eighth_degree: defaultdict= defaultdict(list)

    for key1 in seventh_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(seventh_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(eighth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    eighth_degree[' '.join(sorted(names + [key2]))].append(seventh_degree[key1][0].intersection(set(all_names[key2])))

    
        
    nineth_degree: defaultdict= defaultdict(list)

    for key1 in eighth_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(eighth_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(nineth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    nineth_degree[' '.join(sorted(names + [key2]))].append(eighth_degree[key1][0].intersection(set(all_names[key2])))

    
    tenth_degree: defaultdict= defaultdict(list)

    for key1 in nineth_degree.keys():
        for key2 in all_names.keys():
            if key2 not in key1:
                names = key1.split(' ')
                if (list(nineth_degree[key1][0].intersection(set(all_names[key2]))) != [] and len(tenth_degree[' '.join(sorted(names + [key2]))]) == 0):
                    tenth_degree[' '.join(sorted(names + [key2]))].append(nineth_degree[key1][0].intersection(set(all_names[key2])))
    
    print('\n\n tenth_degree')
    colours = []
    nodesss = []
    nodes_soft = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])
    quantity = 0
    for keys, it in tenth_degree.items():
        for key in keys.split(' '):
            
            if key in c1:
                colours.append(key)
                nodesss.append('green')
            elif key in c2:
                colours.append(key)
                nodesss.append('red')
            elif key in c3:
                colours.append(key)
                nodesss.append('blue')
            else:
                colours.append(key)
                nodesss.append('yellow')
            nodes_soft.loc[quantity] = [key, str(list(it[0])[0]), 'Undirected', 1]
            quantity += 1
        if str(list(it[0])[0]) in c1:
            colours.append(str(list(it[0])[0]))
            nodesss.append('green')
        elif str(list(it[0])[0]) in c2:
            colours.append(str(list(it[0])[0]))
            nodesss.append('red')
        elif str(list(it[0])[0]) in c3:
            colours.append(str(list(it[0])[0]))
            nodesss.append('blue')
        else:
            colours.append(str(list(it[0])[0]))
            nodesss.append('yellow')
        print(f"{keys}: {str(list(it[0])[0])}")
    

    fig = matplotlib.pyplot.figure(figsize=(12,12))
    G = nx.from_pandas_edgelist(nodes_soft, edge_attr=True)
    nx.draw_kamada_kawai(G, ax=fig.add_subplot(), font_size=8, with_labels=True, nodelist=colours, node_size=200, edge_color='orange', node_color=nodesss)
    matplotlib.use("Agg")
    fig.savefig('tenth' + '.png')

    colours = []
    nodesss = []
    nodes_soft = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])
    quantity = 0
    for keys, it in nineth_degree.items():
        for key in keys.split(' '):
            if key not in nodesss:
                if key in c1:
                    colours.append(key)
                    nodesss.append('green')
                elif key in c2:
                    colours.append(key)
                    nodesss.append('red')
                elif key in c3:
                    colours.append(key)
                    nodesss.append('blue')
                else:
                    colours.append(key)
                    nodesss.append('yellow')
            for i in list(it[0]):
                nodes_soft.loc[quantity] = [key, i, 'Undirected', 1]
            quantity += 1
    

    fig = matplotlib.pyplot.figure(figsize=(12,12))
    G = nx.from_pandas_edgelist(nodes_soft, edge_attr=True)
    nx.draw_kamada_kawai(G, ax=fig.add_subplot(), font_size=8, with_labels=True, nodelist=colours, node_size=200, edge_color='orange', node_color=nodesss)
    matplotlib.use("Agg")
    fig.savefig('nineth' + '.png')
    
    
    ################## all rna
     
    rez_dict: defaultdict= defaultdict(int)

    nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])
    quantity = 0
    separator = ''
    not_found = []
    for gene in gene_list:
        entered = 0
        for file in fantom_files:
            
            if gene in file:
                entered = 1
                if len(file) > 3:
                    i=0
                    graph_file = ''
                    for junc in file:
                        
                        if 0 != i % 2 or len(file) == i:
                            graph_file += junc + ","
                        else: 
                            graph_file += junc

                else:
                    graph_file = f"{args.phantom_folder}/{separator.join(file)}"
                graph_results = pd.read_csv(graph_file, delimiter=",", skiprows=1)

                for i in range(len(graph_results)):
                    if (graph_results.iloc[i, 3] >= .95 and not graph_results.iloc[i].isnull().iloc[-3]):
                        rez_dict[' '.join(sorted([gene, graph_results.iloc[i, -3]]))] += 1
                    else:
                        break
        if entered == 0:
            not_found.append(gene)

    nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])

    for key, it in rez_dict.items():
        wan, two = key.split(' ')
        nodes.loc[quantity] = [wan, two, 'Undirected', it]
        quantity += 1

    fig = matplotlib.pyplot.figure(figsize=(40,40))
    G = nx.from_pandas_edgelist(nodes, source='source', target='target', edge_attr='weight')
    nx.draw_kamada_kawai(G, ax=fig.add_subplot(), font_size=8, with_labels=True, node_size=100, edge_color='orange')
    matplotlib.use("Agg")
    fig.savefig('all_rna_graph' + '.png')
     
    ################## all protein
    
    rez_dict: defaultdict= defaultdict(int)

    nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])
    quantity = 0
    separator = ''
    not_found = []
    for gene in prot_list:
        entered = 0
        for file in fantom_files:
            
            if gene in file:
                entered = 1
                if len(file) > 3:
                    i=0
                    graph_file = ''
                    for junc in file:
                        
                        if 0 != i % 2 or len(file) == i:
                            graph_file += junc + ","
                        else: 
                            graph_file += junc

                else:
                    graph_file = f"{args.phantom_folder}/{separator.join(file)}"
                graph_results = pd.read_csv(graph_file, delimiter=",", skiprows=1)

                for i in range(len(graph_results)):
                    if (graph_results.iloc[i, 3] >= .95 and not graph_results.iloc[i].isnull().iloc[-3]):
                        rez_dict[' '.join(sorted([gene, graph_results.iloc[i, -3]]))] += 1
                    else:
                        break
        if entered == 0:
            not_found.append(gene)

    nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])

    for key, it in rez_dict.items():
        wan, two = key.split(' ')
        nodes.loc[quantity] = [wan, two, 'Undirected', it]
        quantity += 1

    fig = matplotlib.pyplot.figure(figsize=(40,40))
    G = nx.from_pandas_edgelist(nodes, source='source', target='target', edge_attr='weight')
    nx.draw_kamada_kawai(G, ax=fig.add_subplot(), font_size=8, with_labels=True, node_size=100, edge_color='orange')
    matplotlib.use("Agg")
    fig.savefig('all_prot_graph' + '.png')
    

        # if path == paths_df.iloc[j, 0].split(':')[1]:
        #     nodes = pd.DataFrame(columns=['source', 'target', 'type', 'weight'])
        #     quantity = 0
        #     for file in fantom_files:
        #         for j in range(len(paths_df)):

        #             if file[1] == paths_df.iloc[j, 1]:
                        
                        
                        

        #                 k = 0
        #                 uniques = []
                    

                    

        #                         quantity += 1
        #                     k+=1

        #     if not nodes.empty:
        #         if len(files_path_gene_dict[path]) > 1:
        #             print(path)
        #         fig = matplotlib.pyplot.figure()
        #         G = nx.from_pandas_edgelist(nodes, source='source', target='target', edge_attr='weight')
        #         nx.draw_random(G, ax=fig.add_subplot(), font_size=5, with_labels=True, node_size=150, edge_color='orange')
        #         matplotlib.use("Agg")
        #         fig.savefig('GO:' + path + '.png')
                    
                    
if __name__  == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--path_file', metavar='path_file', required=True, help='path file')

    parser.add_argument('--phantom_folder', metavar='phantom_folder', required=True, help='phantom folder')
    
    parser.add_argument('--results_file', metavar='results_file', required=True, help='results file')
    
    args = parser.parse_args()
    
    fantom_files = []
    
    for (dirpath, dirnames, filenames) in walk(args.phantom_folder):
        for filename in filenames:
            list_rez: list = []

            if ".csv" not in filename or "p1@" not in filename or 'MUTED' in filename:
                continue

            if "," in filename:
                
                filename_split = filename.split(sep=",")
                for gene in filename_split:
                
                    gene_split = gene.split(sep="@")
                
                    list_rez.append(gene_split[0] + "@")
                    list_rez.append(gene_split[1].split(".")[0])
                
                list_rez.append("." + filename_split[-1].split(".")[1])
            else:
                filename_split = filename.split(sep="@")
                    
                list_rez.append(filename_split[0] + "@")
                list_rez.append(filename_split[1].split(".")[0])
                list_rez.append("." + filename_split[1].split(".")[1])
            fantom_files.extend([list_rez])
        
    main(args, fantom_files)