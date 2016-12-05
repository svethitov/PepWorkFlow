'''Ploting functions for PepWorkFlow'''

from random import randint
from plotly.offline import plot
import plotly.graph_objs as go
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from pepwork.clustering import listtodict

def plot3dscatter(recsdict, filename, xaxis=1, yaxis=2, zaxis=3):
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
        xvalues = \
            [item[1][xaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]
        yvalues = \
            [item[1][yaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]
        zvalues = \
            [item[1][zaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]
        color = 'rgb(' + str(randint(0, 255)) + ',' + str(randint(0, 255)) + \
                    ',' + str(randint(0, 255)) + ')'
        textvalues = [item[0] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]

        if cluster < 0:
            size = 3
        else:
            size = 4
            mesh = go.Mesh3d(
                x=xvalues,
                y=yvalues,
                z=zvalues,
                color=color,
                opacity=0.05,
                alphahull=2,
                name='Cluster: {}'.format(cluster),
            )
            data.append(mesh)

        trace = go.Scatter3d(
            x=xvalues,
            y=yvalues,
            z=zvalues,
            mode='markers',
            marker=dict(
                color=color,
                size=size
            ),
            text=textvalues,
            name='Cluster: {}'.format(cluster),
            opacity=0.8
        )
        data.append(trace)

    axisnames = {
        1: 'Peptide length',
        2: 'Helix',
        3: 'Strand',
        4: 'Turn',
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
    plot(fig, filename=filename)

def hist(vector, filename):
    '''Plots histogram using plotly'''
    data = [go.Histogram(x=vector, autobinx=False, xbins=dict(start=min(vector),
                                                              size=5,
                                                              end=(max(vector) + 25)),
                         marker=dict(color='ADD8E6'))]
    fig = go.Figure(data=data)
    plot(fig, filename=filename)

def clustersdendrogram(data, labels, nodecolor):
    '''Plots dendrogram of the hierarchical clustering'''
    def color_func(index, nodecolor=nodecolor):
        if index in nodecolor.keys():
            return nodecolor[index]
        else:
            return '#D3D3D3'
    plt.figure()
    plt.title('Dendrogram')
    plt.xlabel('Samples')
    plt.ylabel('Distance')
    dendrogram(
        data,
        labels=labels,
        link_color_func=color_func
    )
    plt.show()
