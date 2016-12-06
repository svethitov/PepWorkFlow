'''Ploting functions for PepWorkFlow'''

from plotly.offline import plot
import plotly.graph_objs as go
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram

def _get_trace_mesh(cluster, idx, xaxis, yaxis, zaxis):
    '''Makes mesh and 3dscatter for plot3dscatter'''
    if idx == -1:
        name = 'Outliers'
        size = 2
    else:
        name = 'Cluster: {}'.format(idx + 1)
        size = 4

    xvalues = [cluster.ssfeatures[key][xaxis - 1] for key in cluster.records.keys()]
    yvalues = [cluster.ssfeatures[key][yaxis - 1] for key in cluster.records.keys()]
    zvalues = [cluster.ssfeatures[key][zaxis - 1] for key in cluster.records.keys()]
    name = name
    mesh = go.Mesh3d(
        x=xvalues,
        y=yvalues,
        z=zvalues,
        color=cluster.color,
        opacity=0.05,
        alphahull=2,
        name=name
    )
    trace = go.Scatter3d(
        x=xvalues,
        y=yvalues,
        z=zvalues,
        mode='markers',
        marker=dict(
            color=cluster.color,
            size=size
        ),
        text=[name for name in cluster.records.keys()],
        name=name,
        opacity=0.8
    )
    return mesh, trace

def plot3dscatter(clusters, filename, outliers=None, xaxis=1, yaxis=2, zaxis=3):
    '''Given dictionary or list filled with features vectors plot a 3d scatter plot3dscatterplot
    You should also provide features for the different axes.
    Mapping for axes:
    1 - Peptide Length
    2 - Helix
    3 - Strand
    4 - Turn
    5 - Disulfid bridges #
    Uses the plotly library'''

    data = []

    for idx, cluster in enumerate(clusters):
        mesh, trace = _get_trace_mesh(cluster, idx, xaxis, yaxis, zaxis)
        data.append(mesh)
        data.append(trace)

    if outliers is not None:
        mesh, trace = _get_trace_mesh(outliers, -1, xaxis, yaxis, zaxis)
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
        title='3D Plot of Features',
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
        '''Callable function for node color'''
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
