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

def treetraversal_clustering(node, records: OrderedDict, sim_threshold: float=20.0):
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

    print('Creating new tcoffee directory in {} ...'.format(os.getcwd()))
    os.mkdir('tcoffee')
    os.chdir('tcoffee')

    writefasta(currentrecords, 'working.fasta')
    #command = 't_coffee -seq working.fasta -mode accurate -pdb_type dn \
    #           -pdb_min_sim 50 -pdb_min_cov 20'

    # Runs default t-coffee alignment
    command = 't_coffee -seq working.fasta -outfile default62.aln -output clustalw'
    subprocess.run(command.split())

    # Runs default t-coffee alignment with blosum40
    command = 't_coffee -seq working.fasta -matrix blosum40mt \
                -outfile default40.aln -output clustalw'
    subprocess.run(command.split())

    # Runs PSI-coffee alignment
    command = 't_coffee -seq working.fasta -mode psicoffee -outfile psi62.aln -output clustalw'
    subprocess.run(command.split())

    # Runs PSI-coffee alignment with blosum40
    command = 't_coffee -seq working.fasta -matrix blosum40mt -mode psicoffee \
                -outfile psi40.aln -output clustalw'
    subprocess.run(command.split())

    # Runs M-Coffee
    command = 't_coffee -seq working.fasta -mode mcoffee -outfile m62.aln -output clustalw'
    subprocess.run(command.split())

    # Runs M-Coffee with blosum40
    command = 't_coffee -seq working.fasta -matrix blosum40mt -mode mcoffee \
                -outfile m40.aln -output clustalw'
    subprocess.run(command.split())

    # Runs Local alignment using T-Coffee
    command = 't_coffee -seq working.fasta -method lalign_id_pair \
              -outfile local62.aln -output clustalw'
    subprocess.run(command.split())

    # Runs Local alignment using T-Coffee with blosum40
    command = 't_coffee -seq working.fasta -matrix blosum40mt -method \
              lalign_id_pair -outfile local40.aln -output clustalw'
    subprocess.run(command.split())

    ## Runs Expresso alignment
    #command = 't_coffee -seq working.fasta -mode accurate \
    # -pdb_type dn -outfile accurate62.aln -output clustalw'
    #subprocess.run(command.split())

    ## Runs Expresso alignment with blosum40
    #command = 't_coffee -seq working.fasta -matrix blosum40mt -mode accurate -pdb_type dn\
    #  -outfile accurate40.aln -output clustalw'
    #subprocess.run(command.split())

    # Comparing alignments
    with open('alignments.txt', 'w') as myfile:
        myfile.write('./default62.aln\n')
        myfile.write('./psi62.aln\n')
        myfile.write('./m62.aln\n')
        myfile.write('./local62.aln\n')
        #myfile.write('./accurate62.aln\n')
        myfile.write('./default40.aln\n')
        myfile.write('./psi40.aln\n')
        myfile.write('./m40.aln\n')
        myfile.write('./local40.aln\n')
        #myfile.write('./accurate40.aln\n')

    command = 'trimal -compareset alignments.txt -out best.aln'
    subprocess.run(command.split())

    # Trimming the best alignment

    # Getting sequence similarities
    command = 't_coffee -other_pg seq_reformat -in best.aln -output sim'
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
        result = treetraversal_clustering(
            node=node.left,
            records=records,
            sim_threshold=sim_threshold)
        if result is not None:
            clusters += result
        result = treetraversal_clustering(
            node=node.right,
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

        ## Making accurate alignment using Expresso for the cluster
        #command = 't_coffee -seq working.fasta -mode psicoffee'
        #subprocess.run(command.split())

        # Reads the MSA from the file and stores it in a object
        with open('best.aln', 'r') as msafile:
            msa = AlignIO.read(msafile, 'clustal', alphabet=generic_protein)

        # Uses TrimAI Strict procedure to make trimmed msa and read
        command = 'trimal -in best.aln -out trimmed.aln -htmlout html.html -strict'
        subprocess.run(command.split())
        with open('trimmed.aln', 'r') as trimmed:
            trimmed_msa = AlignIO.read(trimmed, 'clustal', alphabet=generic_protein)

        clusters.append(Cluster(thisclusterrecords, _gen_random_color(), msa, trimmed_msa,
                                node.id, idxs))

    # Deletes the temporaly folder tree
    os.chdir(os.pardir)
    shutil.rmtree('./tcoffee')

    return clusters

def _motif_builder(node, records, e_val, joint=False):
    '''Builds motifs using MEME'''
    print('Working in node {} ...'.format(node.id))

    idxs = node.pre_order()
    keyslist = [key for key in records.keys()]
    name = keyslist[idxs[0]]
    print('Number of children {}'.format(len(idxs)))
    currentrecords = []
    for idx in idxs:
        currentrecords.append(records[keyslist[idx]])

    writefasta(currentrecords, 'working.fasta')

    # Purges redundant sequences
    command = 't_coffee -other_pg seq_reformat -in working.fasta \
              -action +trim _seq_%%80'
    print('Removing redundant sequences ...')
    with open('trimmed.fasta', 'w') as trimmedfile:
        subprocess.run(command.split(), stdout=trimmedfile)

    # Search from motifs using MEME
    if not joint:
        command = 'meme trimmed.fasta -mod oops -nmotifs 10 -maxw 100 \
                    -evt {} -maxiter 1000 -o meme_{}_{}'.format(\
                    e_val, name, node.id)
    else:
        command = 'meme trimmed.fasta -mod oops -nmotifs 10 -maxw 100 \
                    -evt {} -maxiter 1000 -o meme_{}_{}'.format(\
                    e_val, node.left, node.right.id)
    subprocess.run(command.split())


def treetraversal_motifs(node, records: OrderedDict, parentnodes: list, e_val: str):
    '''Recursively traverses the tree and finds motifs along the way'''
    print()
    print('Node number: {}'.format(node.id))

    if node.is_leaf():
        print('Node is a leaf ...')
        return

    print('Parent of {} and {}'.format(node.left.id, node.right.id))
    print('Under this node are {} sequences ...'.format(len(node.pre_order())))

    if node.left.id in parentnodes:
        print()
        print('Reached Cluster on node {} ...'.format(node.left.id))
        if len(node.left.pre_order()) > 2:
            _motif_builder(node.left, records, e_val)
        else:
            print('Cluster is composed just of 2 sequences! Continue ...')
    else:
        treetraversal_motifs(node.left, records, parentnodes, e_val)

    if node.right.id in parentnodes:
        print()
        print('Reached Cluster on node {} ...'.format(node.right.id))
        if len(node.right.pre_order()) > 2:
            _motif_builder(node.right, records, e_val)
        else:
            print('Cluster is composed just of 2 sequences! Continue ...')
    else:
        treetraversal_motifs(node.right, records, parentnodes, e_val)

    if node.left.id in parentnodes and node.right.id in parentnodes:
        print()
        print('Reached node with two children clusters ...')
        _motif_builder(node, records, e_val, joint=True)
