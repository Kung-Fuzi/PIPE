# -*- coding: utf-8 -*-

"""
Created on Wed Mar 23 11:54:24 2022
Created by @Kung-Fuzi

PDB pre-preformatting step for PIPE.
This script takes as input a PDB file and outputs a PIPE-formatted PDB file. 
Please make sure the script is on the same path as your PDB file.

Usage:
    python pdb_preformat.py <pdb file>

Example:
    python pdb_preformat.py 10329_paratope.pdb
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
            
        f = args[0] # Input file
        
    # No input file
    else:
        sys.stderr.write(__doc__)
        sys.exit(1)
        
    return f


def run(fhandle):
    """Return a PIPE-formatted PDB file"""
    
    # Delete all lines except ATOM/HETATM
    records = ('ATOM','HETATM','TER','END')
    
    for line in fhandle:
        if line.startswith(records):
            if line.startswith('ENDMDL'):
                continue
            else:
                prev_line = line
                yield line
    
    # Make TER and END lines if necessary
    if not prev_line.startswith('TER'):
        prev_line = 'TER\n'
        yield prev_line
    if not prev_line.startswith('END'):
        prev_line = 'END\n'
        yield prev_line


def main():
    # Check input
    pdbfn = check_input(sys.argv[1:])
    
    # Open file
    pdbfh = open(pdbfn,'r')
    
    # Run epitope pre-processing step on input PDB
    newfh = run(pdbfh)
    
    # Write list to file
    newfn = os.path.splitext(pdbfn)[0]
    
    with open(f'{newfn}_format.pdb','w') as newfile:
        for line in newfh:
            newfile.write(line)
    
    # Close file
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()