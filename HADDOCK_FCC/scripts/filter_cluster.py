# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 11:05:34 2022
Created by @Kung-Fuzi

Cluter filtering step for PIPE.
This script takes as input predicted epitope/paratope text files and a HADDOCK-FCC 
cluster file and output a filtered cluster file

Usage:
    python3.8 filter_cluster.py <epitope1residues.txt> <epitope2residues.txt> 
    <epitope3residues.txt> <paratoperesidues.txt> <namedclusters.txt>

Example:
    python3.8 filter_cluster.py LILRB1_epi1_residues.txt LILRB1_epi2_residues.txt 
    LILRB1_epi3_residues.txt 10329_paratope_residues.txt namedclusters.txt
"""

import os
import sys
import Bio
from Bio.PDB import PDBParser


# Functions
def check_input(args):
    """"Checks input from command line arguments"""
    
    # Reading from input files
    if len(args) == 5:
        for arg in args:
            if not os.path.isfile(arg):
                sys.stderr.write(f'File not found or not readable: {args[0]}\n')
                sys.stderr.write(__doc__)
                sys.exit(1)
        
        epi1fn = args[0]
        epi2fn = args[1]
        epi3fn = args[2]
        parafn = args[3]
        clusfn = args[4]
        
    # No input file
    else:
        sys.stderr.write(__doc__)
        sys.exit(1)
        
    return epi1fn,epi2fn,epi3fn,parafn,clusfn

def parse_tope(topefile):
    """Parse residues from epitope/paratope text file"""
    
    tope = []
    with open(topefile,'r') as f:
        for fline in f:
            residue = fline.rstrip('\n')
            tope.append(residue)

    return tope

def parse_clusters(clusterfile):
    """Parse clusters from """
    
    clusters = []
    with open(clusterfile,'r') as f:
        flines = f.readlines()

    for fline in flines:
        cluster = fline.strip().split()
        del cluster[0]
        del cluster[1]
        clusters.append(cluster)

    return clusters

def filter_clusters(clusters,cutoff,epitopes,paratope,fcc):
    """Filter clusters by epitope/paratope residue restraints"""
    
    filteredclusters = []
    
    # Find epitope/paratope for cluster representatives
    for cluster in clusters:
    
        # Create PDBParser object
        parser = PDBParser()
        
        # Return structure of PDB file
        clusterrep = f'{cluster[1]}.pdb'
        structure = parser.get_structure('structure',clusterrep)
        
        clusterepitope = []
        clusterparatope = []
        clusterepitopefcc = []
        clusterparatopefcc = []
        
        # Find contact residues
        for residue1 in structure[0]['A'].get_residues():
            for residue2 in structure[0]['B'].get_residues():
                try:
                    if (residue1['CA'] - residue2['CA']) < cutoff:
                        if residue1 not in clusterparatope:
                            resname = residue1.get_resname()
                            reschain = residue1.get_full_id()[2]
                            resseq = residue1.get_full_id()[3][1]
                            clusterparatope.append(f'{resname}.{reschain}.{resseq}')
                        if residue2 not in clusterepitope:
                            resname = residue2.get_resname()
                            reschain = residue2.get_full_id()[2]
                            resseq = residue2.get_full_id()[3][1]
                            clusterepitope.append(f'{resname}.{reschain}.{resseq}')
                except KeyError:
                    continue
        
        # Calculate epitope FCC
        for epitope in epitopes:
            predictedepitope = set(epitope)
            commoncontacts = [res for res in clusterepitope if res in predictedepitope]
            clusterepitopefcc.append(len(commoncontacts)/len(epitope))
    
        # Calculate paratope FCC
        predictedparatope = set(paratope)
        commoncontacts = [res for res in clusterparatope if res in predictedparatope]
        clusterparatopefcc.append(len(commoncontacts)/len(paratope))
        
        # Filter clusters based on FCC
        if clusterparatopefcc[0] > fcc:
            if clusterepitopefcc[0] > fcc or clusterepitopefcc[1] > fcc or clusterepitopefcc[2] > fcc:
                filteredclusters.append(cluster)
    
    return filteredclusters

def main():
    # Check input
    epi1fn,epi2fn,epi3fn,parafn,clusfn = check_input(sys.argv[1:])
    
    # Run cluster filtering step
    epi1 = parse_tope(epi1fn)
    epi2 = parse_tope(epi2fn)
    epi3 = parse_tope(epi3fn)
    epis = [epi1,epi2,epi3]
    para = parse_tope(parafn)
    clus = parse_clusters(clusfn)
    filteredclusters = filter_clusters(clus,7.5,epis,para,0.25)
    
    with open('filteredclusters.txt','a') as f:
        for cluster in filteredclusters:
            line = ' '.join(cluster) + '\n'
            f.write(line)


if __name__ == '__main__':
    main()