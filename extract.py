'''Extract AC from the .list file
connects to UniProt downloads the record and
saves all the records in Biopyhton library'''


import statistics
import pickle
from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt


def getrecords(filename):
    '''Extract the list of AC and returns them in a set of SwissProt Record objects'''
    print('Read list...')
    aclist = []
    with open(filename, 'r') as myfile:
        for line in myfile:
            aclist.append(line.rstrip('\n'))
    print('Completed')
    print('Fetching records...')
    records = set()
    for record in aclist:
        print('Fetching record {} ...'.format(record))
        handle = ExPASy.get_sprot_raw(record)
        records.add(SwissProt.read(handle))
    print('Completed')
    return records

def getpdb(recordlist):
    '''Extract all records with pdb cross-reference'''
    pdbrecords = set()
    for record in recordlist:
        for crossref in record.cross_references:
            if crossref[0] == 'PDB':
                print('PDB cross-reference(s) for ID: {} AC: {} :'.\
                    format(record.entry_name, record.accessions))
                print(crossref)
                pdbrecords.add(record)
    print('{} records found to have PDB cross-reference'.format(len(pdbrecords)))
    return pdbrecords

def getseqstat(seqrecords):
    '''Gives basic statistics for the sequences in iterable object
    Uses matplotlyb.pyplot'''
    print('Getting sequences lengths...')
    lens = [len(seq.sequence) for seq in seqrecords]
    print('Basic stats:')
    print('Mean: {}'.format(statistics.mean(lens)))
    print('Standard deviation: {}'.format(statistics.pstdev(lens)))
    print('Maximum length: {}'.format(max(lens)))
    print('Minimum length: {}'.format(min(lens)))
    print('Range: {}'.format(max(lens)-min(lens)))
    print(sorted(lens))
    plt.hist(lens, bins='auto')
    plt.title('Histogram of Sequences lengths')
    plt.show()

def gettrimmedseq(seqrecords, minlength, maxlength):
    '''Returns set of SwissProt Records filtered by sequence length.
    Range is inclusive'''
    return {record for record in seqrecords if len(record.sequence) >= minlength \
        and len(record.sequence) <= maxlength}

def savebinary(filename, objecttosave):
    '''Saves binary version of the object using pickle highest protocol'''
    with open(filename, 'wb') as myfile:
        pickle.dump(objecttosave, myfile, pickle.HIGHEST_PROTOCOL)

def loadbinary(filename):
    '''Loads binary version of the object using pickle'''
    with open(filename, 'rb') as myfile:
        return pickle.load(myfile)

def writefasta(seqrecords, filename):
    '''Writes records to fasta file ready for alignment'''
    workingseqs = []
    for record in seqrecords:
        print('Preparing {} for writing to file...')
        workingseqs.append(SeqRecord(Seq(record.sequence), id=record.accessions[0]))
    print('Writing fasta file ...')
    SeqIO.write(workingseqs, filename, 'fasta')
    print('Completed')

def settodict(seqset):
    '''Converts set to dict'''
    seqdict = {}
    for record in seqset:
        seqdict[record.accessions[0]] = record
    return seqdict

def subsetseqs(pdbseqsset, clusterdict):
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

    pdbseqdict = settodict(pdbseqsset)

    tofile = []
    for idx, cluster in enumerate(clusters):
        tofile.append([])
        for record in cluster:
            tofile[idx].append(pdbseqdict[record])

    for idx, listfile in enumerate(tofile):
        writefasta(listfile, 'cluster' + str(idx) + '.fasta')
