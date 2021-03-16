import os, sys, gzip, argparse, re, glob

from Bio.PDB import *
from Bio import pairwise2

def read_pdb_files(pdb_files, options_verbose):
    print("hola")
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
            for residue in chain:
                if residue in chain:
                    if residue.id[0] != ' ':
                        removeable.append(residue.id)

        #Now that heteroatoms are selected, remove them from the chain
            for residue in removeable:
                chain.detach_child(residue)
            chains.append(chain)

        #Finally, obtain the alpha carbon chain and store it
        #Check if the length of the polypeptide chain is long enough to not be considered a ligand/cofactor
            if len(chain)<=25:
                structure[0].detach_child(chain.id)
            else:
                chain_alpha = alpha_carbons.build_peptides(chain)
                alpha_carbon_chains.append(chain_alpha[0].get_sequence())


        #Check the file only contains two chains(pairwise interaction) after all the parsing:
        if len(alpha_carbon_chains)!= 2:
            if options_verbose:
                sys.stderr.write("File %s is not a pairwise interaction" % (file))
            continue

        chain_identity = sequence_alignment(alpha_carbon_chains[0], alpha_carbon_chains[1])

        #check binary interactions and add to homo- or heterodimer dictionary (95% similarity threshold)
        if chain_identity > 0.95:
            homodimer_dict[id]=structure
        else:
            heterodimer_dict[id]=structure

    if (len(homodimer_dict)==0) and (len(heterodimer_dict)==0):
        sys.stderr.write("No binary interactions were found")
        exit()
    elif (len(homodimer_dict)==0):
        dict_to_return["heterodimers"]=heterodimer_dict
    elif (len(heterodimer_dict)==0):
        dict_to_return["homodimers"]=homodimer_dict
    else:
        dict_to_return["heterodimers"]=heterodimer_dict
        dict_to_return["homodimers"]=homodimer_dict
    return dict_to_return

def sequence_alignment(chain1,chain2):
    """Comparing if the pairwise interaction holds a homodimer or heterodimer"""
    align=pairwise2.align.globalxx(chain1,chain2)
    identity=align[0][2]/max(len(chain1),len(chain2))
    return identity


    alignment = pairwise2.align.globalxx(sequence1, sequence2)
    return alignment

def dir_path(string):
    """A function to check whether a string is a directory or not"""
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def transform_to_structure(dictionary,name):
    """
    From a dictionary of binary interactions (heterodimer/homodimer dictionaries)
    transform it to a structure class object which will be used in superimposition"""

    structure_object=Structure.Structure(name)

    i=0

    for structure_chains in dictionary.values():
        structure_object.add(Model.Model(i))
        for chain in structure_chains.get_chains():
            structure_object[i].add(chain)
        i+=1

    return structure_object

def check_files(path):
    """
    A function to check whether PDB input files have correct format
    """
    work_files=[]
    my_pattern=re.compile("\w+_\w+_\w+.pdb*")
    for file in os.listdir(path):
        work_files.append(my_pattern.match(file).group(0))
    if not work_files:
#    if my_pattern.match(file) == None:
        raise ValueError("Check the PDB input files format")
    else:
        os.chdir(path)
        return work_files

def output_dir(string, options_force):
    """
    A function to check whether outputfile already exists
    """
    if  options_force is False:
        if os.path.isdir(string):
            raise ValueError("Directory already exists. Please set -f to True to overwrite the directory")
        else:
            sys.stderr.write("Setting the output directory to %s" % (string))
            os.mkdir(string)

def superimpose_structure(fixed_structure, moving_structure, RMSD):
    """
    A function to superimpose two structures using Superimposer from Biopython
    INPUT: fixed structure as refernce, moving structure as the one to be imposed
    """

    sup = Superimposer()   # Superimposer from Biopython
    fixed_atoms = {}
    moving_atoms = {}
    for chainy in fixed_structure:
        for chainy2 in moving_structure:
            fixed_atoms = chain.get_atoms()     # get list of atoms from fixed reference structure
            moving_atoms = chain2.get_atoms()   # get list of atoms from non-fixed structure
            sup.set_atoms(fixed_atoms, moving_atoms) # using list of atoms to create rotation and translation matrix


    rotation = sup.rotran # get rotation from Superimposer
    RMSD = sup.rms # get RMSD from superimposer
    print(sup.rotran)
    print(sup.rms)

    superimposed_structure = sup.apply(moving) # apply the rotation and translation matrix to the non-fixed structure

# # Specify the atom lists
# # 'fixed' and 'moving' are lists of Atom objects
# # The moving atoms will be put on the fixed atoms
# sup.set_atoms(fixed, moving)
# # Print rotation/translation/rmsd
# print(sup.rotran)
# print(sup.rms)
# # Apply rotation/translation to the moving atoms
# sup.apply(moving)
