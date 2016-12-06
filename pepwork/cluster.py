'''Cluster class'''

from collections import OrderedDict
import pepwork.extract

class Cluster:
    '''Class holding the records in a single cluster'''
    def __init__(self, records, color, parentnode=None, membernodes=None):
        self.records = records
        self.color = color
        self.parentnode = parentnode
        self.membernodes = membernodes
        self.ssfeatures = OrderedDict()

    def save_fasta(self, filename: str):
        '''Writes fasta file with the records'''
        pepwork.extract.writefasta(self.records, filename=filename)

    def get_stats(self, filename: str):
        '''Displays basic sequence stats'''
        pepwork.extract.getseqstat(self.records, filename=filename)

    def get_keys(self) -> list:
        '''Returns a list with keys for this cluster'''
        return [key for key in self.records.keys()]
