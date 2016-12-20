'''Extract AC from the .list file
connects to UniProt downloads the record and
saves all the records in Biopyhton library'''

import statistics
import pickle
import time
import urllib
import copy
from collections import OrderedDict
from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pepwork.plots import hist
from pepwork.clustering import toint

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
        while True:
            try:
                time.sleep(0.1)
                print('Fetching record {} ...'.format(record))
                handle = ExPASy.get_sprot_raw(record)
                records.add(SwissProt.read(handle))
                break
            except (urllib.error.HTTPError, urllib.error.URLError):
                print('Network error wating for a second ...')
                time.sleep(1)

    print('All records fetched ...')
    return itertodict(records)

def _get_start_stop(feature_start, feature_stop):
    '''Catches ValueErrors if position unknown'''
    try:
        start_aa = toint(feature_start)
    except ValueError:
        print('There was an error in reading the begining AA.')
        start_aa = int(input('Please enter the value manually: '))

    try:
        end_aa = toint(feature_stop)
    except ValueError:
        print('There was an error in reading the end AA.')
        end_aa = int(input('Please enter the value manually: '))

    return start_aa, end_aa


def trim_sequence(records: OrderedDict) -> OrderedDict:
    '''Returns the same dictionary with SwissProt Objects but with
    trimmed sequences only to the CHAIN or PEPTIDE region'''
    print('Creating trimmed version of the sequences ...')
    output_dict = copy.deepcopy(records)

    for key in output_dict.keys():
        start_aa = None
        print('Finding begining and end of sequence {} ...'.format(key))
        for feature in output_dict[key].features:
            if feature[0] == 'CHAIN':
                print('CHAIN annotation found ...')
                start_aa, end_aa = _get_start_stop(feature[1], feature[2])
                print('Begins at {} ...'.format(start_aa))
                print('Ends at {} ...'.format(end_aa))
                break
        if start_aa is None:
            for feature in output_dict[key].features:
                if feature[0] == 'PEPTIDE':
                    print('PEPTIDE annotation found ...')
                    start_aa, end_aa = _get_start_stop(feature[1], feature[2])
                    print('Begins at {} ...'.format(start_aa))
                    print('Ends at {} ...'.format(end_aa))
                    break

        output_dict[key].sequence = output_dict[key].sequence[start_aa - 1:end_aa]

    return output_dict


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

def writefasta(seqrecords: OrderedDict, filename):
    '''Writes records to fasta file ready for alignment'''
    workingseqs = []
    for record in seqrecords.keys():
        print('Preparing {} for writing to file...'.format(seqrecords[record].accessions[0]))
        workingseqs.append(SeqRecord(Seq(seqrecords[record].sequence),
                                     id=seqrecords[record].accessions[0])
                          )
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
