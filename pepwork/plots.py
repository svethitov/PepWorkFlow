'''Ploting functions for PepWorkFlow'''

from plotly.offline import plot
import plotly.graph_objs as go
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram

def _get_trace(cluster, idx):
    if idx == -1:
        name = 'Outliers'
        size = 2
    else:
        name = 'Cluster: {}'.format(idx + 1)
        size = 4

    return go.Scatter3d(
        x=cluster.ssfeatures_df.Lenght,
        y=cluster.ssfeatures_df.Alpha,
        z=cluster.ssfeatures_df.Beta,
        mode='markers',
        marker=dict(
            color=cluster.color,
            size=size
        ),
        text=cluster.ssfeatures_df.index,
        name=name,
        opacity=0.8
    )

def _get_x(clusters, outliers, idx) -> list:
    '''Get the x axis values'''
    data = []
    for cluster in clusters:
        data.append(cluster.ssfeatures_df.loc[:, idx])

    if outliers is not None:
        data.append(outliers.ssfeatures_df.loc[:, idx])

    return data

def plot3dscatter(clusters, filename, outliers=None):
    '''Given dictionary or list filled with features vectors plot a 3d scatter plot3dscatterplot
    You should also provide features for the different axes.
    Uses the plotly library'''
    data = []

    for idx, cluster in enumerate(clusters):
        data.append(_get_trace(cluster, idx))

    if outliers is not None:
        data.append(_get_trace(outliers, -1))

    layout = go.Layout(
        margin=dict(
            l=1,
            r=1,
            b=1,
            t=1
        ),
        title='3D Plot of Features',
        #scene=dict(
        #    xaxis=dict(
        #        title=clusters[0].ssfeatures_df.columns[0]
        #    ),
        #    yaxis=dict(
        #        title=clusters[0].ssfeatures_df.columns[1]
        #    ),
        #    zaxis=dict(
        #        title=clusters[0].ssfeatures_df.columns[2]
        #    )
        #),
        updatemenus=list([
            dict(
                x=0,
                y=0.85,
                buttons=list([
                    dict(
                        args=['x', _get_x(clusters, outliers, 'Lenght')],
                        label='x: Lenght',
                        method='restyle'
                    ),
                    dict(
                        args=['x', _get_x(clusters, outliers, 'Alpha')],
                        label='x: Alpha',
                        method='restyle'
                    ),
                    dict(
                        args=['x', _get_x(clusters, outliers, 'Beta')],
                        label='x: Beta',
                        method='restyle'
                    ),
                    dict(
                        args=['x', _get_x(clusters, outliers, 'Turn')],
                        label='x: Turn',
                        method='restyle'
                    ),
                    dict(
                        args=['x', _get_x(clusters, outliers, 'Disulfid')],
                        label='x: Disulfid',
                        method='restyle'
                    )
                ])
            ),
            dict(
                x=0,
                y=0.9,
                buttons=list([
                    dict(
                        args=['y', _get_x(clusters, outliers, 'Lenght')],
                        label='y: Lenght',
                        method='restyle'
                    ),
                    dict(
                        args=['y', _get_x(clusters, outliers, 'Alpha')],
                        label='y: Alpha',
                        method='restyle'
                    ),
                    dict(
                        args=['y', _get_x(clusters, outliers, 'Beta')],
                        label='y: Beta',
                        method='restyle'
                    ),
                    dict(
                        args=['y', _get_x(clusters, outliers, 'Turn')],
                        label='y: Turn',
                        method='restyle'
                    ),
                    dict(
                        args=['y', _get_x(clusters, outliers, 'Disulfid')],
                        label='y: Disulfid',
                        method='restyle'
                    )
                ])
            ),
            dict(
                x=0,
                y=0.95,
                buttons=list([
                    dict(
                        args=['z', _get_x(clusters, outliers, 'Lenght')],
                        label='z: Lenght',
                        method='restyle'
                    ),
                    dict(
                        args=['z', _get_x(clusters, outliers, 'Alpha')],
                        label='z: Alpha',
                        method='restyle'
                    ),
                    dict(
                        args=['z', _get_x(clusters, outliers, 'Beta')],
                        label='z: Beta',
                        method='restyle'
                    ),
                    dict(
                        args=['z', _get_x(clusters, outliers, 'Turn')],
                        label='z: Turn',
                        method='restyle'
                    ),
                    dict(
                        args=['z', _get_x(clusters, outliers, 'Disulfid')],
                        label='z: Disulfid',
                        method='restyle'
                    )
                ])
            )
        ])
    )
    fig = go.Figure(data=data, layout=layout)
    plot(fig, filename)

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
