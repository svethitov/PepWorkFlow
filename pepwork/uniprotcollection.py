'''Defines the classes necessary for the Workflow'''

import os
import shutil
import subprocess
from collections import OrderedDict
import pandas as pd
from pepwork.cluster import Cluster
import pepwork.extract
import pepwork.tree
from pepwork.clustering import getfeaturesvector
from pepwork.plots import clustersdendrogram, plot3dscatter

class UniProtCollection:
    '''Class for storing collection of SwissProt records and manipulating them'''
    def __init__(self, mode='list'):
        '''Constructor for UniprotCollection object'''
        if mode == 'list':
            self.allrecords = pepwork.extract.getrecords()
        elif mode == 'bin':
            self.allrecords = pepwork.extract.loadbinary()

        self.allrecords_trimmed = pepwork.extract.trim_sequence(self.allrecords)
        self.records = pepwork.extract.getpdb(self.allrecords_trimmed)
        self.ssfeatures = getfeaturesvector(self.records)
        self.features_df = pd.DataFrame.from_dict(
            self.ssfeatures, orient='index')
        self.features_df.columns = ['Lenght', 'Alpha', 'Beta', 'Turn', 'Disulfid']
        self.clusters = []
        self.outliers = None
        self.clustersdict = {key: -1 for key in self.records.keys()}
        self.rootnode = None
        self.linkagematrix = None
        self.nodecolor = dict()
        self.parentnodes = []

    def save_records(self):
        '''Saves all records of the object into binary
        file. This file can be used to construct the object'''
        pepwork.extract.savebinary(self.allrecords)

    def get_stats(self, scope='all'):
        '''Print statistics of the records'''
        if scope == 'pdb':
            pepwork.extract.getseqstat(self.records, filename='pdb_records.html')
        else:
            pepwork.extract.getseqstat(self.allrecords, filename='all_records.html')

    def buildguidetree(self, method='average'):
        '''Builds guide tree for clustering'''
        self.rootnode, self.linkagematrix = pepwork.tree.buildtree(self.ssfeatures, method)

    def cluster(self, sim_threshold=20.0):
        '''Clusters the pdb cross-refed records'''
        self.clusters = pepwork.tree.treetraversal_clustering(
            node=self.rootnode,
            records=self.records,
            sim_threshold=sim_threshold)

        for idx, cluster in enumerate(self.clusters):
            for key in cluster.get_keys():
                self.clustersdict[key] = idx

        # Set color for each node dependent on the color of the cluster
        for cluster in self.clusters:
            self.nodecolor[cluster.parentnode] = cluster.color
            for node_idx in cluster.membernodes:
                self.nodecolor[node_idx] = cluster.color

        # Gets all outliers in seperate cluster
        records_keys = set([key for key in self.records.keys()])
        records_in_clusters = set()
        for cluster in self.clusters:
            for key in cluster.records.keys():
                records_in_clusters.add(key)
        outliers = records_keys.difference(records_in_clusters)
        outliers_dict = {key: self.records[key] for key in outliers}
        self.outliers = Cluster(outliers_dict, '#D3D3D3')

        # Gets ss features for each record in cluster
        # from the master ss features dictionary
        # and fill dataframe
        for cluster in self.clusters + [self.outliers]:
            for key in cluster.records.keys():
                cluster.ssfeatures[key] = self.ssfeatures[key]
            cluster.ssfeatures_df = self.features_df\
                .loc[[key for key in cluster.records.keys()]]

        self.outliers.ssfeatures_df = self.features_df\
                .loc[[key for key in outliers]]

        # Gets all nodes that are root for a cluster
        for cluster in self.clusters:
            self.parentnodes.append(cluster.parentnode)

    def find_motifs(self, e_val='0.005'):
        '''Find Motifs for the clusters'''
        ## Create BLAST DB from all sequences without pdb cross-ref

        # Deletes blast directory if present
        if os.path.isdir('./blast'):
            print('Deleting old blast directory ...')
            shutil.rmtree('./blast')

        print('Creating blast directory ...')
        os.mkdir('blast')
        os.chdir('blast')
        os.mkdir('blastdb')
        os.chdir('blastdb')
        keys_without_3d = \
            set([key for key in self.allrecords_trimmed]) - set([key for key in self.records])
        records_without_3d = OrderedDict()
        for key in keys_without_3d:
            records_without_3d[key] = self.allrecords_trimmed[key]
        pepwork.extract.writefasta(records_without_3d, 'no_3d.fasta')

        # Makes blast db
        print('Creating BLAST DataBase ...')
        command = 'makeblastdb -in no_3d.fasta -parse_seqids -dbtype prot -out working'
        subprocess.run(command.split())
        os.environ['BLASTDB'] = os.getcwd()
        os.chdir(os.pardir)

        # Extends each cluster
        for idx, cluster in enumerate(self.clusters):
            cluster.save_msa('Cluster_{}.fasta'.format(idx))
            command = 'psiblast -in_msa Cluster_{}.fasta -ignore_msa_master \
                        -out cluster_{}.csv -outfmt=10 -evalue 0.01 -db working'.format(idx, idx)
            print('Running PSIBLAST with Cluster {}'.format(idx))
            subprocess.run(command.split())
            with open('cluster_{}.csv'.format(idx)) as csv_file:
                cluster.extra_records = pd.read_csv(csv_file, header=None, index_col=1)

        os.chdir(os.pardir)


        os.mkdir('MEME_motifs')
        os.chdir('MEME_motifs')

        for idx, cluster in enumerate(self.clusters):
            print('Finding Motifs for Cluster {}'.format(idx))
            working_records = OrderedDict()
            for key in cluster.get_keys():
                working_records[key] = self.allrecords_trimmed[key]
            for key in cluster.extra_records.index:
                working_records[key] = self.allrecords_trimmed[key]

            pepwork.extract.writefasta(working_records, 'working.fasta')

            # Purges redundant sequences
            command = 't_coffee -other_pg seq_reformat -in working.fasta \
                    -action +trim _seq_%%90'
            print('Removing redundant sequences ...')
            with open('trimmed.fasta', 'w') as trimmedfile:
                subprocess.run(command.split(), stdout=trimmedfile)

            command = 'meme trimmed.fasta -mod zoops -nmotifs 8 -evt {}\
                      -maxiter 1000 -o meme_cluster_{}'.format(e_val, idx)
            subprocess.run(command.split())

        #pepwork.tree.treetraversal_motifs(
        #    node=self.rootnode,
        #    records=self.records,
        #    parentnodes=self.parentnodes,
        #    e_val=e_val
        #)
        os.chdir(os.pardir)


    def plot_dendrogram(self):
        '''Plots dendrogram corresponding to the clustering'''
        if self.linkagematrix is None:
            print('Clustering is not done on this object yet ...')
            return 1

        clustersdendrogram(self.linkagematrix, [key for key in self.ssfeatures.keys()],
                           self.nodecolor)

    def plot_3dscatter(self, filename='3dscatter.html'):
        '''Plots 3D scatterplot for all pdb records'''
        plot3dscatter(
            clusters=self.clusters,
            filename=filename,
            outliers=self.outliers
        )
