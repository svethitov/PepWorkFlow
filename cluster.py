'''Clustering and utils for clustering used for the peptides'''


import numpy as np
from sklearn.cluster import DBSCAN, MeanShift, estimate_bandwidth
from sklearn.preprocessing import robust_scale

def getfeaturesvector(seqrecords):
    '''Extract feature vector from iterable object of SwissProt Records'''
    seqdict = dict()
    for record in seqrecords:
        print('Getting features for record {}'.format(record.accessions[0]))
        peptidelength = 0
        print('Searching for PEPTIDE annotation...')
        for feature in record.features:
            if feature[0] == 'PEPTIDE':
                peptidelength = int(feature[2]) - int(feature[1]) + 1
                print('PEPTIDE length found: {} AA'.format(peptidelength))
                break
        if peptidelength == 0:
            print('Searching for CHAIN annotation...')
            for feature in record.features:
                if feature[0] == 'CHAIN':
                    peptidelength = int(feature[2]) - int(feature[1]) + 1
                    print('CHAIN length found: {} AA'.format(peptidelength))
                    break
        helix, turn, strand = 0, 0, 0
        disulfid = 0
        print('Getting secondary structure...')
        for feature in record.features:
            if feature[0] == 'HELIX':
                helix += (int(feature[2]) - int(feature[1]) + 1)
            elif feature[0] == 'TURN':
                turn += (int(feature[2]) - int(feature[1]) + 1)
            elif feature[0] == 'STRAND':
                strand += (int(feature[2]) - int(feature[1]) + 1)
            elif feature[0] == 'DISULFID':
                disulfid += 1
        print('HELIX: {}, STRAND {}, TURN {}, DISULFID {}'.format(\
            helix, strand, turn, disulfid))
        seqdict[record.accessions[0]] = [peptidelength, helix/peptidelength, strand/peptidelength, \
            turn/peptidelength, disulfid]
    return seqdict

def listtodict(recordlist: list):
    '''Returns dictionary object from iterable'''
    return {item[0]: item[1] for item in recordlist}


def scale(recsdict: dict):
    '''Scales the data using Robust Scale'''
    recordslist = []
    for key, value in recsdict.items():
        recordslist.append([key, value])
    print('Creating numpy array...')
    x_raw = np.array([record[1] for record in recordslist])
    return robust_scale(x_raw), recordslist

def labelclusters(recordslist: list, labels: list):
    '''Returns labeled dictionary with the records'''
    for idx, record in enumerate(recordslist):
        record[1].append(labels[idx])
    return listtodict(recordslist)

def removeprevlabel(recsdict):
    '''Removes previous feature label in the dictionary'''
    print('Removing previous labels, if any ...')
    for value in recsdict.values():
        if len(value) > 5:
            value = value[0:4]

    return recsdict

def clustermeanshift(recsdict: dict):
    '''Clusters proteins by MeanShift Algorithm'''
    recsdict = removeprevlabel(recsdict)
    x_scaled, recordslist = scale(recsdict)
    print('Calculating Bandwidth...')
    bandwidth = estimate_bandwidth(x_scaled, n_jobs=-1)
    means = MeanShift(bandwidth=bandwidth, bin_seeding=True, n_jobs=-1)
    print('Clustering using MeanShift ...')
    means.fit(x_scaled)
    print('Calculating basic stats...')
    labels = means.labels_

    uniquelabels = np.unique(labels)
    nclusters = len(uniquelabels)
    print('Number of estimated clusters: {}'.format(nclusters))

    return labelclusters(recordslist, labels), nclusters

def clusterdbscan(recsdict, eps=0.5, min_samples=5):
    '''Clusters proteins by DBSCAN Algorithm'''
    recsdict = removeprevlabel(recsdict)
    x_scaled, recordslist = scale(recsdict)
    print('Clustering using DBSCAN ...')
    mydbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(x_scaled)
    labels = mydbscan.labels_
    nclusters = len(set(labels))
    print('Number of estimated clusters: {}'.format(nclusters))

    return labelclusters(recordslist, labels)
