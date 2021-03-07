import os, sys, gzip, argparse, re, glob

from Bio.PDB import *
from Bio import pairwise2

def read_pdb_files(pdb_files):
    """Given a pdb file, read it, remove the heteroatoms and create a dictionary with the chain ids and the structure

    Input:
        PDB File (files argument) with a pairwise interaction

    Output:
    Dictionary with three elements: Chain ids (2) and the structure """

    dict_to_return={}
    homodimer_dict={}
    heterodimer_dict={}
    pdb_parser=PDBParser(PERMISSIVE=1, QUIET=True)
    alpha_carbons=CaPPBuilder()

    for file in pdb_files:
        id=file[:-4]
        structure=pdb_parser.get_structure(id,file)
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
        #Check if the length of the polypeptide chain is long enough to not be considered a ligand/cofactor
            if len(chain)<=25:
                structure[0].detach_child(chain)
            else:
                chain_alpha = alpha_carbons.build_peptides(chain)
                alpha_carbon_chains.append(chain_alpha[0].get_sequence())

        #Check the file only contains two chains(pairwise interaction) after all the parsing:
        if len(alpha_carbon_chains)!= 2:
            if options.verbose:
                sys.stderr.write("File %s is not a pairwise interaction" % (file))
            continue

        chain_identity = sequence_alignment(alpha_carbon_chains[0], alpha_carbon_chains[1])

        #check binary interactions and add to homo- or heterodimer dictionary (95% similarity threshold)
        if chain_identity > 0.95:
            homodimer_dict[id]=structure
        else:
            heterodimer_dict[id]=structure

    if isEmpty(homodimer_dict) and isEmpty(heterodimer_dict):
        sys.stderr.write("No binary interaction were found")
        exit()
    elif isEmpty(homodimer_dict):
        return heterodimer_dict
    elif isEmpty(heterodimer_dict):
        return homodimer_dict
    else:
        return homodimer_dict,heterodimer_dict

def sequence_alignment(chain1,chain2):
    """Comparing if the pairwise interaction holds a homodimer or heterodimer"""
    align=pairwise2.align.globalxx(chain1,chain2)
    identity=align[0][2]/max(len(chain1),len(chain2))
    return identity

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
        work_files.append(my_pattern.match(file).group(0))
    if not work_files:
#    if my_pattern.match(file) == None:
        raise ValueError("Check the input files format")
    else:
        os.chdir(path)
        return work_files

def output_dir(string):
    """A function to check whether outputfile already exists"""
    if  options.force is False:
        if os.path.isdir(string):
            raise ValueError("Directory already exists. Please set -f to True to overwrite the directory")
        else:
            sys.stderr.write("Setting the output directory to %s" % (string))
            os.mkdir(string)
