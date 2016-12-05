'''Extract AC from the .list file
connects to UniProt downloads the record and
saves all the records in Biopyhton library'''

import statistics
import pickle
import time
from collections import OrderedDict
from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pepwork.plots import hist

def getrecords(filename='list.list'):
    '''Extract the list of AC and returns them in a set of SwissProt Record objects'''
    print('Read list...')
    aclist = []
    with open(filename, 'r') as myfile:
        for line in myfile:
            aclist.append(line.rstrip('\n'))
    print('Completed ...')
    print('Fetching records ...')
    records = set()
    for record in aclist:
        time.sleep(0.25)
        print('Fetching record {} ...'.format(record))
        handle = ExPASy.get_sprot_raw(record)
        records.add(SwissProt.read(handle))
    print('Completed')
    return itertodict(records)

def getpdb(recordsdict):
    '''Extract all records with pdb cross-reference'''
    pdbrecords = set()
    for key in recordsdict.keys():
        for crossref in recordsdict[key].cross_references:
            if crossref[0] == 'PDB':
                print('PDB cross-reference(s) for ID: {} AC: {} :'.\
                    format(recordsdict[key].entry_name, recordsdict[key].accessions))
                print(crossref)
                pdbrecords.add(recordsdict[key])
    print('{} records found to have PDB cross-reference ...'.format(len(pdbrecords)))
    return itertodict(pdbrecords)

def getseqstat(seqrecords, filename=''):
    '''Gives basic statistics for the sequences in iterable object
    Plots histogram of sequence lengths if filename is provided'''
    print('Getting sequences lengths...')
    lens = [len(seqrecords[key].sequence) for key in seqrecords.keys()]
    print('Basic stats:')
    print('Number of sequences: {}'.format(len(lens)))
    print('Mean: {}'.format(statistics.mean(lens)))
    print('Standard deviation: {}'.format(statistics.pstdev(lens)))
    print('Maximum length: {}'.format(max(lens)))
    print('Minimum length: {}'.format(min(lens)))
    print('Range: {}'.format(max(lens)-min(lens)))
    print(sorted(lens))
    if filename != '':
        hist(lens, filename)

def gettrimmedseq(seqrecords, minlength, maxlength):
    '''Returns set of SwissProt Records filtered by sequence length.
    Range is inclusive'''
    return {record for record in seqrecords if len(record.sequence) >= minlength \
        and len(record.sequence) <= maxlength}

def savebinary(objecttosave, filename='records.bin'):
    '''Saves binary version of the object using pickle highest protocol'''
    with open(filename, 'wb') as myfile:
        pickle.dump(objecttosave, myfile, pickle.HIGHEST_PROTOCOL)

def loadbinary(filename='records.bin'):
    '''Loads binary version of the object using pickle'''
    with open(filename, 'rb') as myfile:
        return pickle.load(myfile)

def writefasta(seqrecords, filename):
    '''Writes records to fasta file ready for alignment'''
    workingseqs = []
    for record in seqrecords:
        print('Preparing {} for writing to file...'.format(record.accessions[0]))
        workingseqs.append(SeqRecord(Seq(record.sequence), id=record.accessions[0]))
    print('Writing fasta file ...')
    SeqIO.write(workingseqs, filename, 'fasta')
    print('File written ...')

def itertodict(seqset: iter):
    '''Converts Set to Ordered Dictionary'''
    seqdict = OrderedDict()
    for record in seqset:
        seqdict[record.accessions[0]] = record
    return seqdict

def _dicttoset(seqsdict):
    '''Converts Dictionary to Set'''
    seqsset = set()
    for key in seqsdict.keys():
        seqsset.add(seqsdict[key])
    return seqsset

def _dicttolist(seqsdict):
    '''Converts Dictionary to list'''
    seqslist = []
    for record in seqsdict.values():
        seqslist.append(record)
    return seqslist

def subsetseqs(pdbseqdict, clusterdict):
    '''Subsets the original seqs set based on clusters'''
    uniqueclusters = set([value[5] for value in clusterdict.values()])
    try:
        uniqueclusters.remove(-1)
    except KeyError:
        print('No outliers in the data')

    clusters = []
    for cluster in uniqueclusters:
        clusters.append([key for key, value in clusterdict.items() \
            if value[5] == cluster])

    tofile = []
    for idx, cluster in enumerate(clusters):
        tofile.append([])
        for record in cluster:
            tofile[idx].append(pdbseqdict[record])

    for idx, listfile in enumerate(tofile):
        writefasta(listfile, 'cluster' + str(idx) + '.fasta')
