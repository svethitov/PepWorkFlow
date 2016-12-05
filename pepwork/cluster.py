'''Cluster class'''

import extract

class Cluster:
    '''Class holding the records in a single cluster'''
    def __init__(self, records):
        self.records = records

    def save_fasta(self, filename: str):
        '''Writes fasta file with the records'''
        extract.writefasta(self.records, filename=filename)

    def get_stats(self, filename: str):
        '''Displays basic sequence stats'''
        extract.getseqstat(self.records, filename=filename)

    def get_keys(self) -> list:
        '''Returns a list with keys for this cluster'''
        return [key for key in self.records.keys()]
