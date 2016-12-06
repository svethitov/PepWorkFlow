'''Defines the classes necessary for the Workflow'''

import pepwork.extract
import pepwork.tree
from pepwork.cluster import Cluster
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

        self.records = pepwork.extract.getpdb(self.allrecords)
        self.ssfeatures = getfeaturesvector(self.records)
        self.clusters = []
        self.outliers = None
        self.clustersdict = {key: -1 for key in self.records.keys()}
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
            pepwork.extract.getseqstat(self.records, filename='pdb_records.html')
        else:
            pepwork.extract.getseqstat(self.allrecords, filename='all_records.html')

    def cluster(self):
        '''Clusters the pdb cross-refed records'''
        self.rootnode, self.linkagematrix = pepwork.tree.buildtree(self.ssfeatures)
        self.clusters = pepwork.tree.treetraversal(
            node=self.rootnode,
            records=self.records)
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

        # Gets ss features for each cluster from the master ss features dictionary
        for cluster in self.clusters + [self.outliers]:
            for key in cluster.records.keys():
                cluster.ssfeatures[key] = self.ssfeatures[key]


    def plot_dendrogram(self):
        '''Plots dendrogram corresponding to the clustering'''
        if self.linkagematrix is None:
            print('Clustering is not done on this object yet ...')
            return 1

        clustersdendrogram(self.linkagematrix, [key for key in self.ssfeatures.keys()],
                           self.nodecolor)

    def plot3dscatter(self, filename='3dscatter.html', xaxis=1, yaxis=2, zaxis=3):
        '''Plots 3D scatterplot for all pdb records'''
        plot3dscatter(
            clusters=self.clusters,
            filename=filename,
            outliers=self.outliers,
            xaxis=xaxis,
            yaxis=yaxis,
            zaxis=zaxis
        )
