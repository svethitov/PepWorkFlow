from pepwork.uniprotcollection import UniProtCollection

mycollection = UniProtCollection('bin')
mycollection.cluster()


import copy
from pepwork.plots import plot3dscatter

mytest = {}

for key in mycollection.ssfeatures:
    mytest[key] = copy.deepcopy(mycollection.ssfeatures[key])
    mytest[key].append(mycollection.clustersdict[key])


plot3dscatter(mytest, 'test.html')