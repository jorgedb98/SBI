import os, sys, gzip, argparse, re, glob

from functions import *
from Bio.PDB import *

parser=argparse.ArgumentParser(description="SBI_PYT program")
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-i','--input',
                    required=True,
                    dest="input",
                    action="store",
                    help='file directory. The files that will be sought for must have the format <name>_<chain1>_<chain2>.pdb. Files can be compressed (gz) or not. ',
                    type=dir_path)

parser.add_argument('-s', '--stechiometry',
                    default=None,
                    dest="stec",
                    type=dir_path,
                    action="store",
                    help='Directory containing the appropiate stechiometry file for the molecular complex')

requiredNamed.add_argument('-o','--output',
                    dest="output",
                    required=True,
                    action="store",
                    type=str,
                    help="Directory for the output files. If the directory already exists, check the -f parameter.")

parser.add_argument('-f','--force',
                    dest="force",
                    default=False,
                    action="store_true",
                    help="In case the output directory exists, set the value to true in order to replace it. Otherwise, the output will be stored in a subdirectory.")

parser.add_argument('-v','--verbose',
                    dest="verbose",
                    default=None,
                    action="store_true",
                    help="Set the value to true if you want to console output of the progress")

options=parser.parse_args()

if __name__=="__main__":

    try:
        work_files=check_files(options.input)
    except NotADirectoryError as e:
        sys.stderr.write("Input option does not correspond to an existing directory. Please try again.")
        exit()

    if options.verbose:
        sys.stderr.write("%d files found in $s \n" % (len(work_files),options.input))

    pdb_parser=PDBPaser(PERMISSIVE=1)

#Get all the files as structures in a dictionary with the pdb parser from the Bio package
    structures={}
    for file in work_files:
        structures[file[-7:-4]]=read_pdb_files(file)
