#!/usr/bin/python
'''Workflow for Peptide clustering and alignment'''

import subprocess
import os
import sys
import shutil
from clustering import getfeaturesvector, clusterdbscan
from tree import buildtree
from extract import getseqstat, savebinary, loadbinary, getrecords, getpdb, subsetseqs
from plots import plot3dscatter

def _writingfastas():
    '''Function invoking file writing functionaloty'''
    # write all records with pdb cross-ref
    #writefasta(pdbswissrecords, 'allpdbrecords.fasta')

    # get all fasta files corresponding to clusters
    print('Reading fastas and starting t-coffee alignment ...')
    fastas = []
    for clusterfile in os.listdir():
        if clusterfile.endswith('.fasta'):
            fastas.append(clusterfile)

    if os.path.isdir('./tcoffee'):
        print('Deleting old tcoffee directory ...')
        shutil.rmtree('./tcoffee')
    print('Creating new tcoffee directory ...')
    os.mkdir('tcoffee')
    os.chdir('tcoffee')
    for cluster in fastas:
        current = cluster.split('.')[0]
        print('Working on {}'.format(current))
        os.mkdir(current)
        os.chdir(current)
        filepath = os.path.join(os.pardir, os.pardir, cluster)
        command = 't_coffee -other_pg seq_reformat -in {} -output sim_idscore'\
                    .format(filepath)
        print('Estimating diversity and writing to file ...')
        with open('sim_idscore.txt', 'w') as simfile:
            subprocess.run(command.split(), stdout=simfile)
        command = 't_coffee -seq {} -mode accurate -pdb_type dn'.format(filepath)
        print('Running alignment ...')
        subprocess.run(command.split())
        command = 't_coffee -other_pg seq_reformat -in {}.aln -output sim'.format(current)
        print('Estimating diversity on the alignment and writing to file ...')
        with open('sim.txt', 'w') as simfile:
            subprocess.run(command.split(), stdout=simfile)
        os.chdir(os.pardir)


def _plotting(clustersdict):
    '''Function invoking plotting functionality'''
    axes = [1, 2, 3]
    filename = '3dscatter.html'
    while True:
        print('')
        print('Plotting 3D graph with the following settings:')
        print('Name: {}'.format(filename))
        print('x: {}, y: {}, z: {}'.format(axes[0], axes[1], axes[2]))
        answer = input('Type "no" if you want to change them:')
        if answer == 'no':
            axes[0] = input('Enter value for first axis: ')
            axes[1] = input('Enter value for second axis: ')
            axes[2] = input('Enter value for third axis: ')
            axes = [int(axis) for axis in axes]
            filename = input('Enter a value for filename: ')

        plot3dscatter(clustersdict, filename=filename, xaxis=axes[0],
                      yaxis=axes[1], zaxis=axes[2])
        answer = input('Are you happy with the plot? Type "y" to continue: ')
        if answer == 'y':
            break

def _clustering(featuresvector, pdbswissrecords):
    '''Function invoking dbscan clustering procedure'''
    eps = 1.1
    min_samples = 3
    length_weight = 1
    ss_weight = 1
    while True:
        print('')
        print('Clustering will be performed with the following settings:')
        print('eps: {}, min_samples: {}, length_weight: {}, ss_weight: {}'
              .format(eps, min_samples, length_weight, ss_weight))
        answer = input('Type "no" if you want to change them:')
        if answer == 'no':
            eps = input('Enter new value for eps: ')
            eps = float(eps)
            min_samples = input('Enter new value for min_samples: ')
            min_samples = int(min_samples)
            length_weight = input('Enter new value for length_weight: ')
            length_weight = float(length_weight)
            ss_weight = input('Enter new value for ss_weight: ')
            ss_weight = float(ss_weight)
        print('Clustering ...')
        clustersdict = clusterdbscan(featuresvector, eps=eps,
                                     min_samples=min_samples, length_weight=length_weight,
                                     ss_weight=ss_weight)

        _plotting(clustersdict=clustersdict)

        answer = input('Are you happy with the result of the clustering? Type "y" to continue: ')
        if answer == 'y':
            break

   # subsetseqs(pdbseqsdict=pdbswissrecords, clusterdict=clustersdict)

    _writingfastas()

def _treebuilding(featuresvector):
    '''Function invoking UPGMA tree building procedure'''
    root_node = buildtree(featuresvector=featuresvector)


def main():
    '''Main function for the WorkFlow'''
    print('Searching for list.list or RecordsSet.bin ...')
    if os.path.isfile('./RecordsSet.bin'):
        print('RecordsSet.bin found ...')
        swissrecords = loadbinary('RecordsSet.bin')
    elif os.path.isfile('./list.list'):
        print('list.list found ...')
        swissrecords = getrecords('list.list')
    else:
        print('No usable file found! Exiting')
        sys.exit()

    getseqstat(swissrecords, 'all_records.html')
    input('Press Enter to continue ...')
    savebinary('RecordsSet.bin', swissrecords)
    pdbswissrecords = getpdb(swissrecords)
    getseqstat(pdbswissrecords, 'pdb_cross-refed_records.html')
    input('Press Enter to continue ...')

    featuresvector = getfeaturesvector(pdbswissrecords)

    print()
    print('By default will use the new tree building procedure ...')
    choice = input('Press Enter to continue or type "no" to change to old dbscan: ')

    if choice == 'no':
        _clustering(featuresvector=featuresvector, pdbswissrecords=pdbswissrecords)
    else:
        _treebuilding(featuresvector=featuresvector)


if __name__ == '__main__':
    main()
