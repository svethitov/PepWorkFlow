#!/usr/bin/python
'''Workflow for Peptide clustering and alignment'''

import subprocess
import os
import sys
import shutil
from cluster import getfeaturesvector, clusterdbscan
from extract import getseqstat, savebinary, loadbinary, getrecords, getpdb, subsetseqs
from plots import plot3dscatter


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

    eps = 1.1
    min_samples = 3
    while True:
        print('')
        print('Clustering will be performed with the following settings:')
        print('eps: {}, min_samples: {}'.format(eps, min_samples))
        answer = input('Type "no" if you want to change them:')
        if answer == 'no':
            eps = input('Enter new value for eps: ')
            eps = float(eps)
            min_samples = input('Enter new value for min_samples: ')
            min_samples = int(min_samples)
        print('Clustering ...')
        clustersdict = clusterdbscan(featuresvector, eps=eps, min_samples=min_samples)

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

        answer = input('Are you happy with the result of the clustering? Type "y" to continue: ')
        if answer == 'y':
            break

    subsetseqs(pdbseqsset=pdbswissrecords, clusterdict=clustersdict)

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


if __name__ == '__main__':
    main()
