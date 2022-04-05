# -*- coding: utf-8 -*-

"""
Created on Wed Mar 23 11:54:24 2022
Created by @Kung-Fuzi

PDB pre-processing step for PIPE.
This script takes as input a PDB file and outputs a minimal PDB file.

Usage:
    python3.8 preprocess_pdb.py <pdb file>

Example:
    python3.8 preprocess_pdb.py 10329_paratope.pdb
"""


import os
import sys


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

def run(fhandle):
    """Return a PIPE-formatted PDB file"""
    
    # Delete all lines except ATOM/HETATM
    modelrecords = ('MODEL','ENDMDL')
    structurerecords = ('TER','END')
    atomrecords = ('ATOM','HETATM')
    
    for line in fhandle:
        if line.startswith(modelrecords):
            continue
        elif line.startswith(structurerecords):
            if line.startswith('TER'):
                prev_line = 'TER\n'
                yield prev_line
            elif line.startswith('END'):
                prev_line = 'END\n'
                yield prev_line
        elif line.startswith(atomrecords):
            prev_line = line
            yield line
    
    # Make TER and END lines if necessary
    if not prev_line.startswith(structurerecords):
        prev_line = 'TER\n'
        yield prev_line
        prev_line = 'END\n'
        yield prev_line

def main():
    # Check input
    pdbfn = check_input(sys.argv[1:])
    
    # Run epitope pre-processing step on input PDB
    inputpdb = open(pdbfn,'r')
    pdbfh = run(inputpdb)
    
    # Write to file
    pdbdir = os.path.dirname(pdbfn)
    pdbbase = os.path.basename(pdbfn)
    newfn = os.path.join(pdbdir,f'new_{pdbbase}')
    with open(newfn,'w') as outputpdb:
        for line in pdbfh:
            outputpdb.write(line)
    
    inputpdb.close()
    sys.exit(0)


if __name__ == '__main__':
    main()