import os, sys, gzip, argparse, re, glob

from Bio.PDB import *

def read_pdb_files(pdb_files):
    """Given a pdb file, read it, remove the heteroatoms and create a dictionary with the chain ids and the structure

    Input:
        PDB File (files argument) with a pairwise interaction

    Output:
    Dictionary with three elements: Chain ids (2) and the structure """

    dict_to_return={}

    for file in pdb_files:
        structure=pdb_parser.get_structure(files[:-4],files)

        chains_ids=''.join([chain.id for chain in structure.get_chains()])
        chains=[]
        alpha_carbon_chains=[]

        #Obtain the alpha carbon structure of each chain
        removeable=[]
        for chain in structure.get_chains():
            chains.append(chain)
            for residue in chain:
                if residue in chain:
                    if residue.id[0] != ' ':
                        removeable.append(residue.id)

        #Now that heteroatoms are selected, remove them from the chain
            for residue in removeable:
                chain.detach_child(residue)

        #Finally, obtain the alpha carbon chain and store it
            chain_alpha = alpha_carbons.build_peptides(chain)
            alpha_carbon_chains.append(chain_alpha[0].get_sequence())

    #Check the file only contains two chains(pairwise interaction) after all the parsing:
    if len(alpha_carbon_chains != 2):
        if options.verbose:
            sys.stderr.write("File %s is not a pairwise interaction" % (file))
        continue

def dir_path(string):
    """A function to check whether a string is a directory or not"""
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def check_files(path):
    """A function to check whether inputfiles have correct format"""
    work_files=[]
    my_pattern=re.compile("\w+_\w+_\w+.pdb*")
    for file in os.listdir(path):
        work_files.append(my_pattern.match(file))
    if not work_files:
#    if my_pattern.match(file) == None:
        raise ValueError("Check the input files format")
    else:
        return work_files

def output_dir(string):
    """A function to check whether outputfile already exists"""
    if  options.force is False:
        if os.path.isdir(string):
            raise ValueError("Directory already exists. Please set -f to True to overwrite the directory")
        else:
            sys.stderr.write("Setting the output directory to %s" % (string))
            os.mkdir(string)
