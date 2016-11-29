'''Clustering and utils for clustering used for the peptides'''

import plotly.plotly as py
from plotly.offline import plot
import plotly.graph_objs as go
from random import randint
import numpy as np
from sklearn.cluster import DBSCAN, MeanShift, estimate_bandwidth
from sklearn.preprocessing import robust_scale

def getfeaturevector(seqrecords):
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
    '''Returns dictionary object from list'''
    return {item[0]: item[1] for item in recordlist}

def plot3dscatter(recsdict: dict or list, xaxis=1, yaxis=2, zaxis=3, n_clusters=1):
    '''Given dictionary or list filled with features vectors plot a 3d scatter plot3dscatterplot
    You should also provide features for the different axes.
    Mapping for axes:
    1 - Peptide Length
    2 - Helix %
    3 - Strand %
    4 - Turn %
    5 - Disulfid bridges #
    Uses the plotly library'''
    if isinstance(recsdict, list):
        recsdict = listtodict(recsdict)

    vectors = [item for item in recsdict.items()]
    data = []

    # Gets lenght of Value
    for key in recsdict.keys():
        mykey = key
        break

    if len(recsdict[mykey]) == 6:
        uniqueclusters = set([item[1][5] for item in vectors])
    else:
        uniqueclusters = [1]

    for cluster in uniqueclusters:
        if cluster < 0:
            size = 4
        else:
            size = 8
        trace = go.Scatter3d(
            x=[item[1][xaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster],
            y=[item[1][yaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster],
            z=[item[1][zaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster],
            mode='markers',
            marker=dict(
                color='rgb(' + str(randint(0, 255)) + ',' + str(randint(0, 255)) + \
                    ',' + str(randint(0, 255)) + ')',
                size=size
            ),
            text=[item[0] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster],
            name='Cluster: {}'.format(cluster),
            opacity=0.8
        )
        data.append(trace)

    axisnames = {
        1: 'Peptide length',
        2: 'Helix %',
        3: 'Strand %',
        4: 'Turn %',
        5: 'Disulfid bridges #'
    }
    layout = go.Layout(
        margin=dict(
            l=1,
            r=1,
            b=1,
            t=1
        ),
        title='3d plot of features',
        scene=dict(
            xaxis=dict(
                title=axisnames[xaxis]
            ),
            yaxis=dict(
                title=axisnames[yaxis]
            ),
            zaxis=dict(
                title=axisnames[zaxis]
            )
        )
    )
    fig = go.Figure(data=data, layout=layout)
    plot(fig, filename='test_graph.html')

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

def clustermeanshift(recsdict: dict):
    '''Clusters proteins by MeanShift Algorithm'''
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
    x_scaled, recordslist = scale(recsdict)
    print('Clustering using DBSCAN ...')
    mydbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(x_scaled)
    labels = mydbscan.labels_
    nclusters = len(set(labels))
    print('Number of estimated clusters: {}'.format(nclusters))

    return labelclusters(recordslist, labels), nclusters
