'''Clustering and utils for clustering used for the peptides'''

from collections import OrderedDict
import numpy as np
from sklearn.cluster import DBSCAN, MeanShift, estimate_bandwidth
from sklearn.preprocessing import robust_scale

def toint(number):
    '''Implements int() function but checks if there is special chraracter in the string'''
    if isinstance(number, str):
        return int(''.join(e for e in number if e.isalnum()))
    elif isinstance(number, float):
        return int(number)
    else:
        return number

def secondarystructuretrim(startaa, endaa, featurestartaa, featureendaa):
    '''Trimmes length of Secondary Structure features if they are outside
    of CHAIN or PEPTIDE sequence'''
    if startaa > featureendaa or endaa < featurestartaa:
        #This are the cases where the feature is outside the region of interest
        return 0

    if featurestartaa < startaa:
        start = startaa
    else:
        start = featurestartaa

    if featureendaa > endaa:
        end = endaa
    else:
        end = featureendaa

    return end - start + 1

def getfeaturesvector(records: OrderedDict) -> OrderedDict:
    '''Constructs features vector for every record based on the secondary structure'''
    vectorsdict = OrderedDict()
    for key in records.keys():
        print('Getting features for record {} ...'.format(records[key].accessions[0]))
        peplength = 0
        print('Searching for CHAIN annotation ...')
        for feature in records[key].features:
            if feature[0] == 'CHAIN':
                peplength = toint(feature[2]) - toint(feature[1]) + 1
                print('CHAIN length found: {} AA'.format(peplength))
                startaa = toint(feature[1])
                endaa = toint(feature[2])
                break
        if peplength == 0:
            print('Searching for PEPTIDE annotation ...')
            for feature in records[key].features:
                if feature[0] == 'PEPTIDE':
                    peplength = toint(feature[2]) - toint(feature[1]) + 1
                    print('PEPTIDE length found: {} AA'.format(peplength))
                    startaa = toint(feature[1])
                    endaa = toint(feature[2])
                    break
        helix, turn, strand, disulfid = 0, 0, 0, 0
        print('Getting secondary structure ...')
        for feature in records[key].features:
            if feature[0] == 'HELIX':
                helix += secondarystructuretrim(startaa=startaa, endaa=endaa,
                                                featurestartaa=toint(feature[1]),
                                                featureendaa=toint(feature[2]))
            elif feature[0] == 'TURN':
                turn += secondarystructuretrim(startaa=startaa, endaa=endaa,
                                               featurestartaa=toint(feature[1]),
                                               featureendaa=toint(feature[2]))
            elif feature[0] == 'STRAND':
                strand += secondarystructuretrim(startaa=startaa, endaa=endaa,
                                                 featurestartaa=toint(feature[1]),
                                                 featureendaa=toint(feature[2]))
            elif feature[0] == 'DISULFID' and feature[1] >= startaa and feature[2] <= endaa:
                disulfid += 1
        print('HELIX: {}, STRAND {}, TURN {}, DISULFID {}'.format(\
            helix, strand, turn, disulfid))

        vectorsdict[records[key].accessions[0]] = [peplength, helix/peplength, strand/peplength, \
            turn/peplength, disulfid]

    return vectorsdict


def listtodict(recordlist):
    '''Returns dictionary object from iterable'''
    return {item[0]: item[1] for item in recordlist}


def scale(recsdict):
    '''Scales the data using Robust Scale'''
    recordslist = []
    for key, value in recsdict.items():
        recordslist.append([key, value])
    print('Creating numpy array...')
    x_raw = np.array([record[1] for record in recordslist])
    return robust_scale(x_raw), recordslist

def labelclusters(recordslist, labels):
    '''Returns labeled dictionary with the records'''
    for idx, record in enumerate(recordslist):
        record[1].append(labels[idx])
    return listtodict(recordslist)

def removeprevlabel(recsdict):
    '''Removes previous feature label in the dictionary'''
    print('Removing previous labels, if any ...')
    for key, value in recsdict.items():
        if len(value) > 5:
            recsdict[key] = value[0:5]

    return recsdict

def clustermeanshift(recsdict):
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
    nclusters = len(uniquelabels) - 1
    print('Number of estimated clusters (without outliers): {}'.format(nclusters))

    return labelclusters(recordslist, labels), nclusters

def clusterdbscan(recsdict, eps=0.5, min_samples=5, length_weight=1, ss_weight=1):
    '''Clusters proteins by DBSCAN Algorithm'''
    recsdict = removeprevlabel(recsdict)
    x_scaled, recordslist = scale(recsdict)
    print('Weighting ...')
    x_scaled[:, 0] *= length_weight
    x_scaled[:, 1] *= ss_weight
    x_scaled[:, 2] *= ss_weight
    x_scaled[:, 3] *= ss_weight
    print('Clustering using DBSCAN ...')
    mydbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(x_scaled)
    labels = mydbscan.labels_
    nclusters = len(set(labels)) - 1
    print('Number of estimated clusters (without outliers): {}'.format(nclusters))

    return labelclusters(recordslist, labels)
