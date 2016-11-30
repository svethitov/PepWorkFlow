#!/usr/bin/python
'''Workflow for Peptide clustering and alignment'''

from cluster import *
from extract import *
from plots import plot3dscatter
import subprocess
import os

def main():
    '''Main function for the WorkFlow'''
    print('The script expect list.list file or PDBRecordsSet.bin:')
    swissrecords = getrecords('list.list')
    getseqstat(swissrecords, 'all_records.html')
    input('Press Enter key to continue ...')
    pdbswissrecords = getpdb(swissrecords)
    getseqstat(pdbswissrecords, 'pdb_cross-refed_records.html')
    input('Press Enter key to continue ...')
    savebinary('PDBRecordsSet.bin', pdbswissrecords)
    featuresvector = getfeaturesvector(pdbswissrecords)

    eps = 1.2
    min_samples = 2
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

    os.mkdir('tcoffee')
    os.chdir('tcoffee')
    for cluster in fastas:
        current = cluster.split('.')[0]
        print('Working on {}'.format(current))
        os.mkdir(current)
        os.chdir(current)
        filepath = os.path.join(os.pardir, os.pardir, cluster)
        command = 't_coffee -seq {} -mode accurate -pdb_type dn'.format(filepath)
        subprocess.run(command.split())
        os.chdir(os.pardir)


if __name__ == '__main__':
    main()
