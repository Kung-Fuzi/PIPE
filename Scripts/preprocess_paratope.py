# -*- coding: utf-8 -*-

"""
Created on Thu Mar 17 13:27:04 2022
Created by @Kung-Fuzi

Paratope pre-processing step for PIPE.
This script takes as an input a preprocessed SAbPred i-Patch output PDB file 
and outputs a text file containing the predicted paratope residues as 
resname.reschain.resseq.

Usage:
    python3.8 preprocess_paratope.py <pdb file>

Example:
    python3.8 preprocess_paratope.py 10329_paratope.pdb
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
            
        fn = args[0] # Input file
        
    # No input file
    else:
        sys.stderr.write(__doc__)
        sys.exit(1)
        
    return fn


def run(inputfile):
    """Return a list of paratope residues"""
    
    # Create PDBParser object
    parser = PDBParser()
    
    # Return structure of PDB file
    recstructure = parser.get_structure('receptor',inputfile)
    
    recparatope = []
    recsafe = []
    recblock = []
    
    # Find all residues within specified B-factor cutoff of i-Patch
    for atom in recstructure.get_atoms():
        if atom.get_bfactor() >= 10: # PIPE i-Patch cutoff
            if atom.get_parent() not in recparatope:
                recparatope.append(atom.get_parent())
    
    # Find residues to block
    for residue1 in recstructure.get_residues():
        for residue2 in recparatope:
            if residue1 != residue2:
                try: # Find all residues within 20A of each paratope residue
                    if (residue1['CA'] - residue2['CA']) < 20:
                        if residue1 not in recsafe:
                            recsafe.append(residue1)
                except KeyError:
                    continue
    
    for residue in recstructure.get_residues():
        if residue not in recsafe:
            recblock.append(residue)
    
    return recparatope, recblock


def main():
    # Check input
    recpdb = check_input(sys.argv[1:])
    
    # Run epitope pre-processing step on input PDB
    para, block = run(recpdb)
    
    # Write list to file
    fn = os.path.splitext(recpdb)[0]
    with open(f'{fn}_residues.txt','w') as newfile1:
        for res in para:
            resname = res.get_resname()
            reschain = res.get_full_id()[2]
            resseq = res.get_full_id()[3][1]
            newfile1.write(f'{resname}.{reschain}.{resseq}\n')
    with open(f'{fn}_blockedresidues.txt','w') as newfile2:
        for res in block:
            resname = res.get_resname()
            reschain = res.get_full_id()[2]
            resseq = res.get_full_id()[3][1]
            newfile2.write(f'{resname}.{reschain}.{resseq}\n')
    sys.exit(0)


if __name__ == '__main__':
    main()