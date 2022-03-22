#!/usr/bin/env python

"""
Outputs the names of the cluster members based on a clustering output file
and the original file listing the PDB files

Authors:
           RODRIGUES Joao

Changed output lines 86-98 to write a namedclusters.txt output file with each 
line containing a string representation of a cluster list with the format:
    ["Cluster 1", "decoy_1", "decoy_2", "decoy_3"]
@jason.kong
"""

import os
import sys

USAGE = "python %s <cluster_x.out> <file.nam>" % os.path.basename(sys.argv[0])


def read_clusters(path):
    """Reads clusters from a FCC output file."""
    clusters = []
    cl_file = open(path, 'r')
    for line in cl_file:
        # Cluster 8 -> 193 141 142 144 151 168 171 172 178
        models = list(map(int, line.split()[3:]))
        clusters.append(models)

    return clusters


def read_list(path):
    """
    Reads a list containing one file per line.
    Returns an index of line number - line content
    """
    with open(path, 'r') as fhandle:
        fdata = {}
        for nline, line in enumerate(fhandle):
            if not line.strip():
                continue
            # Remove extension
            fdata[nline+1] = '.'.join(line.strip().split('.')[:-1])

    return fdata


def cross_data(clusters, flist):
    """Matches names in flist to the numbers in clusters."""
    cluster_l = []
    for cl in clusters:
        ncl = [flist[s] for s in cl]
        cluster_l.append(ncl)

    return cluster_l


if __name__ == '__main__':

    if sys.version_info[0:2] < (3, 8):
        cur_version = "%s.%s" % sys.version_info[0:2]
        sys.stderr.write("- Python version not supported (%s). Please use 3.8 or newer.\n" % cur_version)
        sys.exit(1)

    if len(sys.argv[1:]) != 2:
        print(USAGE)
        sys.exit(1)

    cluster_file = os.path.abspath(sys.argv[1])
    pdblist_file = os.path.abspath(sys.argv[2])
    
    try:
        cl_list = read_clusters(cluster_file)
    except IOError:
        sys.stderr.write('Error: file not found (%s)\nAborting..\n' % cluster_file)
        sys.exit(1)

    try:
        pdb_list = read_list(pdblist_file)
    except IOError:
        sys.stderr.write('Error: file not found (%s)\nAborting..\n' % pdblist_file)
        sys.exit(1)

    named_clusters = cross_data(cl_list, pdb_list)
    
    # Output
    nclines = []
    
    for i, nc in enumerate(named_clusters):
        ncline = "Cluster %i -> %s" % (i+1, ' '.join(nc))
        nclines.append(ncline)
        print(ncline)
    
    with open('namedclusters.txt', 'a') as f:
        for ncline in nclines:
            line = ncline + '\n'
            f.write(line)