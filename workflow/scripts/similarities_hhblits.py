from turtle import width
import matplotlib.pyplot as plt
import pandas as pd
import sys

import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.cm as cm
import matplotlib.colors as colors
from itertools import cycle
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.colors as mc
import colorsys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

def create_dict_length_profiles(df_tsv):
    """Create a dict of dicts with length - profiles from a DataFrame .

    Args:
        df_tsv (pandas.DataFrame): Tsv results from hhblits blasttab output in a dataframe

    Returns:
        dict: dictionnary profile: length
    """
    identical = df_tsv.loc[df_tsv['query'] == df_tsv['target'], :]
    identical = identical.set_index('query').qend.to_dict()

    return identical

##########################################################################

def calculate_coverage(target_start, target_end, target_length):
    """Calculate the coverage of the alignment on the target profile.

    Args:
        target_start (numpy.array): Starts of the profiles
        target_end (numpy.array): Ends of the profiles
        target_length (numpy.array): Length of the profiles

    Returns:
        numpy.array: coverage of the alignment on the target
    """
    return (target_end - target_start +1) / target_length

##########################################################################

def create_adjency_matrix(df_tsv, profiles):
    """Create adjency matrix from a TSV file .

    Args:
        df_tsv (pandas.DataFrame): output format of the hhblits blasttab
        profiles (list): iterable of profiles

    Returns:
        pandas.DataFrame: adjavency matrix
    """
    df_tsv=df_tsv[~(df_tsv["query"] == df_tsv["target"])]

    df = pd.pivot_table(
        data=df_tsv, index="query", columns="target", values="pid", fill_value=0
    )
    idx = df.columns.union(profiles)
    df = df.reindex(index=idx, columns=idx, fill_value=0)

    return df

##########################################################################

def get_color_cmap(name, n_colors=6):

    """
    Return discrete colors from a matplotlib palette.
    :param name: Name of the palette. This should be a named matplotlib colormap.
    :type: str
    :param n_colors: Number of discrete colors in the palette.
    :type: int
    :return: List-like object of colors as hexadecimal tuples
    :type: list
    """

    brewer_qual_pals = {"Accent": 8, "Dark2": 8, "Paired": 12,
                        "Pastel1": 9, "Pastel2": 8,
                        "Set1": 9, "Set2": 8, "Set3": 12, 'tab20':20, 'tab20b':20}


    if name == 'tab20' and n_colors > 19:
        second = 'tab20b'
        ncolor2 = n_colors - 19
        n_colors = 19
    else :
        second = False

    cmap = getattr(cm, name)
    
    if name in brewer_qual_pals:
        bins = np.linspace(0, 1, brewer_qual_pals[name])
        if 'tab20' == name :
            len_bins = len(bins)
            bins = [bins[i] for i in range(len_bins) if i != 14][:n_colors]
        else :
            bins = bins[:n_colors]
    else:
        bins = np.linspace(0, 1, n_colors + 2)[1:-1]

    palette = list(map(tuple, cmap(bins)[:, :3]))

    if second :
        cmap = getattr(cm, second)
        bins = np.linspace(0, 1, brewer_qual_pals[second])[:ncolor2]
        palette += list(map(tuple, cmap(bins)[:, :3]))

        pal_cycle = cycle(palette)
        palette = [next(pal_cycle) for _ in range(n_colors+ncolor2)]
    else :
        pal_cycle = cycle(palette)
        palette = [next(pal_cycle) for _ in range(n_colors)]

    return [colors.rgb2hex(rgb) for rgb in palette]

##########################################################################

def contrasting_text_color(hex_str):
    '''
    Input a string without hash sign of RGB hex digits to compute
    complementary contrasting color such as for fonts
    '''

    (r, g, b) = (hex_str[1:3], hex_str[3:5], hex_str[5:])

    luminosity = 1 - (int(r, 16) * 0.299 + int(g, 16) * 0.587 + int(b, 16) * 0.114) / 255

    return '#131516' if luminosity < 0.5 else 'white'

##########################################################################

def adjust_lightness(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

##########################################################################

def visu_graph(identity_df, output, dict_color={}, threshold=0.5) :
    
    """
    Visualisation of the graph using the weight for the color. If weight >0.5
    put the edge red. And write the graph in graphMl format.
    :params adjacency_mtrix: adjacency matrix calculate with get_matrix_interaction_system()
    :type: pandas.DataFrame
    :params output: Name of the graphml file
    :type: str
    :return: Nothing
    
    """ 
    
    # Create the graph from adjacency matrix
    graph = nx.from_numpy_matrix(identity_df.values)
    graph = nx.relabel_nodes(graph, dict(enumerate(identity_df.columns)))
    outdeg = graph.degree()

    if not snakemake.params.singleton:
        to_remove = [n[0] for n in outdeg if outdeg[n[0]] == 0]

        graph.remove_nodes_from(to_remove)

    # Get name of all nodes/genes
    all_gene = list(graph)
    num_gene = len(all_gene)

    # Create the dict of color if none given
    if not dict_color:
        palette = get_color_cmap("tab20", n_colors = num_gene)
        dict_color = {all_gene[i]:palette[i] for i in range(num_gene)}
    
    # Parsing the edges and changing the color edge to red if percentage id >= 50% (or any threshold)
    # Having a list of the edges with to have the width proportionnal to percentage id
    width_edges = []
    for n1, n2, edge_dict in graph.edges.data() :
        width_edges.append(edge_dict['weight']*10)
        if edge_dict['weight'] >= threshold :
            graph.edges[(n1, n2)]["color"] = "#C41E3A"
        else:
            graph.edges[(n1, n2)]["color"] = '#a9a9a9'


    graph = graph.to_undirected()
    # print("Calculating layout...")    

    # Color edges
    edges,edge_colors = zip(*nx.get_edge_attributes(graph,'color').items())   

    # Choose between : dot, neato, fdp, sfdp, twopi, circo
    pos=graphviz_layout(graph, prog=snakemake.params.graph_layout)

    # Put the color of the node
    nx.set_node_attributes(graph, dict_color, "color")    

    # Color nodes
    nodes,nodes_colors = zip(*nx.get_node_attributes(graph,'color').items())   
    
    # Write the graph
    # nx.write_graphml(graph, output)
 
    # Get all the color for the border if colored border = True
    nodes_edges_color = []
    for node_color in nodes_colors:
        if snakemake.config['default_values_plot']['colored_border']:
            nodes_edges_color.append("#2F3D44" if node_color == "#FFFFFF" else adjust_lightness(node_color))
        else :
            nodes_edges_color.append('#131516')

    plt.figure(figsize=(12,10))

    # If you have too much node maybe reduce the node_size to 50
    nx.draw_networkx_nodes(graph, pos, node_color=nodes_colors, node_size=snakemake.params.size_node, edgecolors=nodes_edges_color)

    nx.draw_networkx_labels(graph,pos, font_size=snakemake.params.font_size)

    nx.draw_networkx_edges(graph,pos,edgelist=edges, edge_color=edge_colors, width=width_edges)

    # Create the legend in case needed
    custom_lines = []
    custom_text = []
    
    for gene, color in dict_color.items() :
        custom_text.append(gene)
        custom_lines.append(Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=color, linestyle='none'))
    
    plt.legend(custom_lines, custom_text, bbox_to_anchor=(1.05, 1), loc='upper left', prop={"size":'small'})
    
    #Label drawing as well
    # nx.draw_networkx_labels(graph,pos,font_size=8)

    plt.axis('off')
    plt.tight_layout()
    
    plt.title(f"Graph of the identity between the profiles using evalue {snakemake.params.e_val}\n and a coverage of {snakemake.params.cov}")
    
    if output!=None: plt.savefig(output, dpi=300, bbox_inches='tight')

    plt.close('all')

    return

##########################################################################################

all_tsv_hhsearch = snakemake.input.tsv

tsv_names = ['query', 'target', 'pid', 'alnLen', 'mismatch', 'gapOpen', 'qstart', 'qend', 'tstart', 'tend', 'e_val', 'score']

tsv_dtypes = {'query':"string",
              'target':"string",
              'pid':"float64",
              'alnLen':"int32",
              'mismatch':"int32",
              'gapOpen':"int32",
              'qstart':"int32",
              'qend':"int32",
              'tstart':"int32",
              'tend':"int32",
              'eval':"float64",
              'score':"float64",
              }

all_tsv = []

for tsv in all_tsv_hhsearch:
    all_tsv.append(pd.read_table(tsv, names=tsv_names, dtype=tsv_dtypes))

# Get all the results in one dataframe
all_tsv = pd.concat(all_tsv, ignore_index=True)

# Get all the length of the profiles
dict_prot = create_dict_length_profiles(df_tsv=all_tsv)

# Get all the coverage of the target
all_tsv['coverage'] = calculate_coverage(target_start=all_tsv.tstart.values,
                                        target_end=all_tsv.tend.values, 
                                        target_length=all_tsv.target.map(dict_prot).values,
                                        )

# Selection of the profiles
all_tsv = all_tsv[(all_tsv.e_val <= snakemake.params.e_val) &
                 (all_tsv.coverage >= snakemake.params.cov) &
                 (all_tsv.pid >= snakemake.params.pid)
                ].reset_index(drop=True)

# Calculate the adjacency matrix
adjacency = create_adjency_matrix(df_tsv=all_tsv,
                                    profiles=dict_prot.keys())

if snakemake.params.annotations:
    color_dict = pd.read_table(snakemake.params.annotations, index_col=0).color.to_dict()
else :
    color_dict = {}

for outfile in snakemake.output:
    visu_graph(adjacency, outfile, color_dict, threshold=0.5)
