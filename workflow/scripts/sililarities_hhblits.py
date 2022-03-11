import matplotlib.pyplot as plt
import pandas as pd
import sys, os

import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.cm as cm
import matplotlib.colors as colors
from itertools import cycle
from matplotlib.lines import Line2D
import numpy as np

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
    all_tsv=all_tsv[~(all_tsv["query"] == all_tsv["target"])]

    df = pd.pivot_table(
        data=df_tsv, index="query", columns="target", values="pid", fill_value=0
    )
    idx = df.columns.union(profiles)
    df = df.reindex(index=idx, columns=idx, fill_value=0)

    return df

##########################################################################

all_tsv_hhsearch = snakemake.input.tsv

tsv_names = ['query', 'target', 'pid', 'alnLen', 'mismatch', 'gapOpen', 'qstart', 'qend', 'tstart', 'tend', 'eval', 'score']

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
all_tsv = all_tsv[(all_tsv.eval <= snakemake.params.e_val) &
                 (all_tsv.coverage >= snakemake.params.cov) &
                 (all_tsv.pid >= snakemake.params.pid)
                ].reset_index(drop=True)

# Calculate the adjacency matrix
adjacency = create_adjency_matrix(df_tsv=all_tsv,
                                    profiles=dict_prot.heys())


