# -*- coding: utf-8 -*-

"""
Created on Thu Mar 17 13:27:04 2022
Created by @Kung-Fuzi

Paratope pre-processing step for PIPE.
This script takes as an input a SAbPred i-Patch output PDB file and outputs a 
text file containing the predicted paratope residues as resname.reschain.resseq.

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
    
    recparatope = []
    recsafe = []
    recblkB = []
    recblkC = []
    
    # Find all residues within specified B-factor cutoff of i-Patch
    for atom in recstructure.get_atoms():
        if atom.get_bfactor() >= 10: # PIPE i-Patch cutoff
            resname = atom.get_parent().get_resname()
            reschain = atom.get_parent().get_full_id()[2]
            resseq = atom.get_parent().get_full_id()[3][1]
            residue = f'{resname}.{reschain}.{resseq}'
            if residue not in recparatope:
                recparatope.append(residue)
    
    # Find residues to block
    for residue in recstructure.get_residues():
        resname = residue.get_resname()
        reschain = residue.get_full_id()[2]
        resseq = residue.get_full_id()[3][1]
        residue1 = f'{resname}.{reschain}.{resseq}'
        for residue2 in recparatope:
            if residue1 != residue2:
                try: # Find all residues within 20A of each paratope residue
                    if (residue1['CA'] - residue2['CA']) < 20:
                        if residue1 not in recsafe:
                            recsafe.append(residue1)
                except KeyError:
                    continue
    
    for residue in recstructure.get_residues():
        resname = residue.get_resname()
        reschain = residue.get_full_id()[2]
        resseq = residue.get_full_id()[3][1]
        residueblk = f'{resname}.{reschain}.{resseq}'
        if residueblk not in recsafe:
            if reschain == 'B':
                recblkB.append(residueblk)
            elif reschain == 'C':
                recblkC.append(residueblk)
    
    return recparatope, recblkB, recblkC


def main():
    # Check input
    recpdb = check_input(sys.argv[1:])
    
    # Run epitope pre-processing step on input PDB
    para, blkb, blkc = run(recpdb)
    
    # Write list to file
    fn = os.path.splitext(recpdb)[0]
    with open(f'{fn}_residues.txt','w') as newfile1:
        for res in para:
            newfile1.write(f'{res}\n')
    with open(f'{fn}_blockedB.txt','w') as newfile2:
        for res in blkb:
            newfile2.write(f'{res}\n')
    with open(f'{fn}_blockedC.txt','w') as newfile3:
        for res in blkc:
            newfile3.write(f'{res}\n')
    sys.exit(0)


if __name__ == '__main__':
    main()