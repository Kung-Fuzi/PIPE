# -*- coding: utf-8 -*-

"""
Created on Wed Mar 23 11:54:24 2022
Created by @Kung-Fuzi

PDB pre-processing step for PIPE.
This script takes as input a PDB file and removes non-ATOM/HETATM info.

Usage:
    python3.8 preprocess_pdb.py <pdb file>

Example:
    python3.8 preprocess_pdb.py 10329_paratope.pdb > temp.pdb
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
    
    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(pdbfh):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken pipes
        pass
    
    # Close file and exit
    inputpdb.close()
    sys.exit(0)


if __name__ == '__main__':
    main()