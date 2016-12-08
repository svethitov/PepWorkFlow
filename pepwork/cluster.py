'''Cluster class'''

from collections import OrderedDict
import pepwork.extract
from Bio import AlignIO

class Cluster:
    '''Class holding the records in a single cluster'''
    def __init__(self, records, color, msa=None, trimmed_msa=None,
                 parentnode=None, membernodes=None):
        self.records = records
        self.color = color
        self.parentnode = parentnode
        self.membernodes = membernodes
        self.ssfeatures = OrderedDict()
        self.msa = msa
        self.trimmed_msa = trimmed_msa

    def save_fasta(self, filename: str):
        '''Writes fasta file with the records'''
        pepwork.extract.writefasta(self.records, filename=filename)

    def get_stats(self, filename: str):
        '''Displays basic sequence stats'''
        pepwork.extract.getseqstat(self.records, filename=filename)

    def get_keys(self) -> list:
        '''Returns a list with keys for this cluster'''
        return [key for key in self.records.keys()]

    def save_msa(self, filename: str, outformat: str='fasta'):
        '''Writes MSA file from the full msa of the cluster'''
        AlignIO.write(self.msa, filename, outformat)

    def save_trimmed_msa(self, filename: str, outformat: str='fasta'):
        '''Writes MSA file from the full msa of the cluster'''
        AlignIO.write(self.trimmed_msa, filename, outformat)
