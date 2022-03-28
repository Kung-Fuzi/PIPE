# -*- coding: utf-8 -*-

"""
Created on Thu Mar 17 13:27:04 2022
Created by @Kung-Fuzi

Paratope pre-processing step for PIPE.
This script takes as an input a SAbPred i-Patch output PDB file and outputs a 
text file containing the predicted paratope residues as resname.resseq.

Usage:
    python3.8 paratope_preprocess.py <pdb file>

Example:
    python3.8 paratope_preprocess.py 10329_paratope.pdb
"""


import os
import sys
import Bio
from Bio.PDB import PDBParser


def check_input(args):
    """"Checks input from command line arguments"""
    
    # Reading from input file
    if len(args) == 1:
        if not os.path.isfile(args[0]):
            sys.stderr.write(f'File not found or not readable: {args[0]}\n')
            sys.stderr.write(__doc__)
            sys.exit(1)
            
        f = args[0] # Input file
        
    # No input file
    else:
        sys.stderr.write(__doc__)
        sys.exit(1)
        
    return f


def run(inputfile):
    """Return a list of paratope residues"""
    
    # Create PDBParser object
    parser = PDBParser()
    
    # Return structure of PDB file
    recstructure = parser.get_structure('receptor',inputfile)
    
    ipatch = []
    recparatope = []
    
    # Find all residues within specified B-factor cutoff of i-Patch
    for atom in recstructure.get_atoms():
        if atom.get_bfactor() >= 10: # PIPE i-Patch cutoff
            residue = atom.get_parent()
            if residue not in ipatch:
                ipatch.append(residue)
                
    for residue1 in recstructure.get_residues():
        for residue2 in ipatch:
            if residue1 != residue2:
                try:
                    if (residue1['CA'] - residue2['CA']) < 10: # PIPE expand cutoff
                        resname = residue1.get_resname()
                        resseq = residue1.get_full_id()[3][1]
                        residue = f'{resname}.{resseq}'
                        if residue not in recparatope:
                            recparatope.append(residue)
                except KeyError:
                    continue
                
    # Sort paratope residues
    recparatopeseqs = []
    recparatopesort = []
    
    for residue in recparatope:
        resseq = int(residue.split('.')[1])
        recparatopeseqs.append(resseq)
    
    recparatopeseqscopy = recparatopeseqs.copy()
    recparatopeseqscopy.sort()
    
    for resseq in recparatopeseqscopy:
        recparatopesort.append(recparatope[recparatopeseqs.index(resseq)])
    
    return recparatopesort


def main():
    # Check input
    recpdb = check_input(sys.argv[1:])
    
    # Run epitope pre-processing step on input PDB
    para = run(recpdb)
    
    # Write list to file
    fn = os.path.splitext(recpdb)[0]
    with open(f'{fn}_residues.txt','w') as newfile:
        for res in para:
            newfile.write(f'{res}\n')
    sys.exit(0)


if __name__ == '__main__':
    main()