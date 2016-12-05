'''Cluster class'''

import pepwork.extract

class Cluster:
    '''Class holding the records in a single cluster'''
    def __init__(self, records, color, parentnode, membernodes):
        self.records = records
        self.color = color
        self.parentnode = parentnode
        self.membernodes = membernodes

    def save_fasta(self, filename: str):
        '''Writes fasta file with the records'''
        pepwork.extract.writefasta(self.records, filename=filename)

    def get_stats(self, filename: str):
        '''Displays basic sequence stats'''
        pepwork.extract.getseqstat(self.records, filename=filename)

    def get_keys(self) -> list:
        '''Returns a list with keys for this cluster'''
        return [key for key in self.records.keys()]
