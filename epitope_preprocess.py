# -*- coding: utf-8 -*-

"""
Created on Thu Mar 17 13:27:04 2022
Created by @Kung-Fuzi

Epitope pre-processing step for PIPE.
This script takes as input a SAbPred EpiPred output PDB file and outputs a 
text file containing the predicted epitope residues. Please make sure the script 
is on the same path as your PDB file.

Usage:
    python epitope_preprocess.py <pdb file>

Example:
    python epitope_preprocess.py LILRB1_epi1.pdb
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
    """Return a list of epitope residues"""
    
    # Create PDBParser object
    parser = PDBParser()
    
    # Return structure of PDB file
    ligstructure = parser.get_structure('ligand',inputfile)
    
    ligepitope = []
    
    # Find all residues within specified B-factor cutoff of EpiPred
    for atom in ligstructure.get_atoms():
        if atom.get_bfactor() == 100: # EpiPred sets epitope residue B-factor to 100
            resname = atom.get_parent().get_resname()
            resseq = atom.get_parent().get_full_id()[3][1]
            residue = f'{resname}.{resseq}'
            if residue not in ligepitope:
                ligepitope.append(residue)
                
    return ligepitope


def main():
    # Check input
    ligpdb = check_input(sys.argv[1:])
    
    # Run epitope pre-processing step on input PDB
    epi = run(ligpdb)
    
    # Write list to file
    fn = os.path.splitext(ligpdb)[0]
    with open(f'{fn}_residues.txt','w') as newfile:
        for res in epi:
            newfile.write(f'{res}\n')