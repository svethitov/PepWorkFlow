'''Ploting functions for PepWorkFlow'''

from random import randint
from plotly.offline import plot
import plotly.graph_objs as go
from cluster import listtodict

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
        uniqueclusters = set(1)

    for cluster in uniqueclusters:
        xvalues = \
            [item[1][xaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]
        yvalues = \
            [item[1][yaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]
        zvalues = \
            [item[1][zaxis - 1] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]

        if cluster < 0:
            size = 4
        else:
            size = 8
            mesh = go.Mesh3d(
                x=xvalues,
                y=yvalues,
                z=zvalues,
                color='87CEFA',
                opacity=0.1
            )
            data.append(mesh)

        textvalues = [item[0] for item in vectors if len(item[1]) < 6 or item[1][5] == cluster]
        trace = go.Scatter3d(
            x=xvalues,
            y=yvalues,
            z=zvalues,
            mode='markers',
            marker=dict(
                color='rgb(' + str(randint(0, 255)) + ',' + str(randint(0, 255)) + \
                    ',' + str(randint(0, 255)) + ')',
                size=size
            ),
            text=textvalues,
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
    plot(fig, filename=filename)

def hist(vector, filename):
    '''Plots histogram using plotly'''
    data = [go.Histogram(x=vector, xbins=dict(
        start=0,
        size=25,
        end=max(vector)
    ))]
    fig = go.Figure(data=data)
    plot(fig, filename=filename)
