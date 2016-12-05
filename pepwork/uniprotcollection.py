'''Defines the classes necessary for the Workflow'''

import pepwork.extract
import pepwork.tree
from pepwork.clustering import getfeaturesvector
from pepwork.plots import clustersdendrogram

class UniProtCollection:
    '''Class for storing collection of SwissProt records and manipulating them'''
    def __init__(self, mode='list'):
        '''Constructor for UniprotCollection object'''
        if mode == 'list':
            self.allrecords = pepwork.extract.getrecords()
        elif mode == 'bin':
            self.allrecords = pepwork.extract.loadbinary()

        self.pdbrecords = pepwork.extract.getpdb(self.allrecords)
        self.ssfeatures = getfeaturesvector(self.pdbrecords)
        self.clusters = []
        self.clustersdict = {key: -1 for key in self.pdbrecords.keys()}
        self.rootnode = None
        self.linkagematrix = None
        self.nodecolor = dict()

    def save_records(self):
        '''Saves all records of the object into binary
        file. This file can be used to construct the object'''
        pepwork.extract.savebinary(self.allrecords)

    def get_stats(self, scope='all'):
        '''Print statistics of the records'''
        if scope == 'pdb':
            pepwork.extract.getseqstat(self.pdbrecords, filename='pdb_records.html')
        else:
            pepwork.extract.getseqstat(self.allrecords, filename='all_records.html')

    def cluster(self):
        '''Clusters the pdb cross-refed records'''
        self.rootnode, self.linkagematrix = pepwork.tree.buildtree(self.ssfeatures)
        self.clusters = pepwork.tree.treetraversal(
            node=self.rootnode,
            records=self.pdbrecords)
        for idx, cluster in enumerate(self.clusters):
            for key in cluster.get_keys():
                self.clustersdict[key] = idx

        for cluster in self.clusters:
            self.nodecolor[cluster.parentnode] = cluster.color
            for node_idx in cluster.membernodes:
                self.nodecolor[node_idx] = cluster.color



    def plot_dendrogram(self):
        '''Plots dendrogram corresponding to the clustering'''
        if self.linkagematrix is None:
            print('Clustering is not done on this object yet ...')
            return 1

        clustersdendrogram(self.linkagematrix, [key for key in self.ssfeatures.keys()],
                           self.nodecolor)
