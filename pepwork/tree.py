'''Adds UPGMA tree building'''

import os
import subprocess
import re
import random
import shutil
from collections import OrderedDict
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import cophenet, to_tree
from scipy.spatial.distance import pdist
from Bio import AlignIO
from Bio.Alphabet import generic_protein
from pepwork.clustering import removeprevlabel, scale
from pepwork.extract import writefasta, itertodict
from pepwork.cluster import Cluster

def _gen_random_color():
    '''Generates random HEX color'''
    return '#' + ''.join([random.choice('0123456789ABCDEF') for x in range(6)])


def buildtree(featuresvector, method):
    '''Creates tree from peptide features and returns root node'''
    featuresvector = removeprevlabel(featuresvector)
    x_scaled, _ = scale(featuresvector)
    print('Building linkage matrix using {} algorithm ...'.format(method))
    linkage_matrix = linkage(x_scaled, method)
    coph, _ = cophenet(linkage_matrix, pdist(x_scaled))
    print('Cophenet parameter (values close to 1 are good): {}'.format(coph))
    return to_tree(linkage_matrix), linkage_matrix


def childrentraversal(node):
    '''Returns the idxs of the children nodes'''
    if node.is_leaf():
        return [node.id]
    else:
        idxs = [node.id]
        idxs += childrentraversal(node.left)
        idxs += childrentraversal(node.right)
        return idxs

def treetraversal(node, records: OrderedDict, sim_threshold: float=20.0):
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

    if os.path.isdir('./tcoffee'):
        print('Deleting old tcoffee directory ...')
        shutil.rmtree('./tcoffee')

    print('Creating new tcoffee directory ...')
    os.mkdir('tcoffee')
    os.chdir('tcoffee')

    writefasta(currentrecords, 'working.fasta')
    #command = 't_coffee -seq working.fasta -mode accurate -matrix=blosum50mt -pdb_type dn \
    #           -pdb_min_sim 50 -pdb_min_cov 20'
    command = 't_coffee -seq working.fasta -outfile=working.aln'
    subprocess.run(command.split())

    command = 't_coffee -other_pg seq_reformat -in working.aln -output sim'
    with open('sim.txt', 'w') as simfile:
        subprocess.run(command.split(), stdout=simfile)

    sim_below_threshold = False
    with open('sim.txt', 'r') as simfile:
        for line in simfile:
            if line.startswith('BOT'):
                similarity = re.search(r'\d{1,2}\.\d{2}', line)
                print('Similarity {}'.format(similarity.group(0)))
                if float(similarity.group(0)) < sim_threshold:
                    print('Pair of records with lower than threshold similarity found ...')
                    sim_below_threshold = True
                    break

    if sim_below_threshold:
        result = treetraversal(node=node.left,
                               records=records,
                               sim_threshold=sim_threshold)
        if result is not None:
            clusters += result
        result = treetraversal(node=node.right,
                               records=records,
                               sim_threshold=sim_threshold)
        if result is not None:
            clusters += result
    else:
        print('Cluster found on node {} with {} peptides ...'.\
            format(node.id, len(idxs)))
        idxs = []
        idxs += childrentraversal(node.left)
        idxs += childrentraversal(node.right)
        thisclusterrecords = itertodict(currentrecords)

        # Making accurate alignment using Expresso for the cluster
        command = 't_coffee -seq working.fasta -mode psicoffee'
        subprocess.run(command.split())

        # Reads the MSA from the file and stores it in a object
        with open('working.aln', 'r') as msafile:
            msa = AlignIO.read(msafile, 'clustal', alphabet=generic_protein)

        # Uses TrimAI Strict procedure to make trimmed msa and read
        command = 'trimal -in working.aln -out trimmed.aln -htmlout html.html -strict'
        subprocess.run(command.split())
        with open('trimmed.aln', 'r') as trimmed:
            trimmed_msa = AlignIO.read(trimmed, 'clustal', alphabet=generic_protein)

        clusters.append(Cluster(thisclusterrecords, _gen_random_color(), msa, trimmed_msa,
                                node.id, idxs))

    # Deletes the temporaly folder tree
    os.chdir(os.pardir)
    shutil.rmtree('./tcoffee')

    return clusters
