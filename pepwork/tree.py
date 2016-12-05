'''Adds UPGMA tree building'''

import os
import subprocess
import re
from collections import OrderedDict
from scipy.cluster.hierarchy import average
from scipy.cluster.hierarchy import cophenet, to_tree
from scipy.spatial.distance import pdist
from clustering import removeprevlabel, scale
from extract import writefasta, itertodict
from cluster import Cluster

def buildtree(featuresvector):
    '''Creates tree from peptide features and returns root node'''
    featuresvector = removeprevlabel(featuresvector)
    x_scaled, _ = scale(featuresvector)
    print('Building linkage matrix ...')
    linkage_matrix = average(x_scaled)
    coph, _ = cophenet(linkage_matrix, pdist(x_scaled))
    print('Cophenet parameter (values close to 1 are good): {}'.format(coph))
    input('Press Enter to continue ...')
    return to_tree(linkage_matrix)

def treetraversal(node, records: OrderedDict):
    '''Recursively traversing the tree and searching for clusters of peptides'''
    print()
    print('Node number: {}'.format(node.id))
    if node.is_leaf():
        print('Node is a leaf ...')
        return
    clusters = []

    idxs = node.pre_order()
    keyslist = [key for key in records.keys()]
    print('Number of children {}'.format(len(idxs)))
    currentrecords = []
    for idx in idxs:
        currentrecords.append(records[keyslist[idx]])

    writefasta(currentrecords, 'working.fasta')
    command = 'clustalo -i working.fasta -o clustalo.aln'
    subprocess.run(command.split())

    command = 't_coffee -other_pg seq_reformat -in clustalo.aln -output sim'
    with open('clustalo_sim.txt', 'w') as simfile:
        subprocess.run(command.split(), stdout=simfile)

    sim_below_threshold = False
    with open('clustalo_sim.txt', 'r') as simfile:
        for line in simfile:
            if line.startswith('BOT'):
                similarity = re.search(r'\d{1,2}\.\d{2}', line)
                print('Similarity {}'.format(similarity.group(0)))
                if float(similarity.group(0)) < 20.0:
                    print('Pair of records with lower than threshold similarity found ...')
                    sim_below_threshold = True
                    break
    os.remove('working.fasta')
    os.remove('clustalo.aln')
    os.remove('clustalo_sim.txt')
    if sim_below_threshold:
        result = treetraversal(node=node.left,
                               records=records)
        if result is not None:
            clusters += result
        result = treetraversal(node=node.right,
                               records=records)
        if result is not None:
            clusters += result
    else:
        print('Cluster found on node {} with {} peptides ...'.\
            format(node.id, len(idxs)))
        thisclusterrecords = itertodict(currentrecords)
        clusters.append(Cluster(thisclusterrecords))
    
    return clusters
