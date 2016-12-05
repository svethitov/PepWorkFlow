'''Defines the classes necessary for the Workflow'''

import extract
import tree
from clustering import getfeaturesvector

class UniProtCollection:
    '''Class for storing collection of SwissProt records and manipulating them'''
    def __init__(self, mode='list'):
        '''Constructor for UniprotCollection object'''
        if mode == 'list':
            self.allrecords = extract.getrecords()
        elif mode == 'bin':
            self.allrecords = extract.loadbinary()

        self.pdbrecords = extract.getpdb(self.allrecords)
        self.ssfeatures = getfeaturesvector(self.pdbrecords)
        self.clusters = []
        self.clustersdict = {key: -1 for key in self.pdbrecords.keys()}
        self.rootnode = None

    def save_records(self):
        '''Saves all records of the object into binary
        file. This file can be used to construct the object'''
        extract.savebinary(self.allrecords)

    def get_stats(self, scope='all'):
        '''Print statistics of the records'''
        if scope == 'pdb':
            extract.getseqstat(self.pdbrecords, filename='pdb_records.html')
        else:
            extract.getseqstat(self.allrecords, filename='all_records.html')

    def cluster(self):
        '''Clusters the pdb cross-refed records'''
        self.rootnode = tree.buildtree(self.ssfeatures)
        self.clusters = tree.treetraversal(node=self.rootnode,
                                           records=self.pdbrecords)
        for idx, cluster in enumerate(self.clusters):
            for key in cluster.get_keys():
                self.clustersdict[key] = idx

