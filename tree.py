'''Adds UPGMA tree building'''

from scipy.cluster.hierarchy import average
from scipy.cluster.hierarchy import cophenet, to_tree
from scipy.spatial.distance import pdist
from cluster import removeprevlabel, scale

def buildtree(featuresvector):
    '''Creates tree from peptide features and returns root node'''
    featuresvector = removeprevlabel(featuresvector)
    x_scaled, _ = scale(featuresvector)
    print('Building linkage matrix ...')
    linkage_matrix = average(x_scaled)
    coph, _ = cophenet(linkage_matrix, pdist(x_scaled))
    print('Cophenet parameter (values close to 1 are good): {}'.format(coph))
    return to_tree(linkage_matrix)
