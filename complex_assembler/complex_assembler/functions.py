##########################################
##  FUNCTIONS SCRIPT FOR MAIN PROGRAM   ##
##########################################

import os, sys, gzip, argparse, glob, string, shutil

from Bio.PDB import *
from Bio import pairwise2



def read_pdb_files(pdb_files, options_verbose):
    """
    Given a pdb file, read it, remove the heteroatoms and create a dictionary with the chain ids and the structure. This dictionary also classifies the
    input in protein-protein interaction, DNA-protein interaction or RNA-protein interaction (this last one, considering if RNA is single-strand or double-strand).

    - Input: PDB File (files argument) with a pairwise interaction

    - Output: Dictionary with three elements: Chain ids (2) and the structure
    """
    # Create empty variables

    dict_with_sNP = {}
    dict_with_dNP = {}
    dict_with_PP = {}

    pdb_parser=PDBParser(PERMISSIVE=1, QUIET=True)

    for file in pdb_files:                                          # Iterate over the files
        id=file[:-4]                                                # Id will be the first 4 characters of the file name
        structure=pdb_parser.get_structure(id,file)                 # Get structure
        chains_ids=''.join([chain.id for chain in structure.get_chains()])
        chains=[]
        alpha_carbon_chains=0                                       # Initialize counter to know if alpha carbons have been added

        # Obtain the alpha carbon structure of each chain
        removeable=[]
        for chain in structure.get_chains():
            for residue in chain:
                if residue in chain:
                    if residue.id[0] != ' ':
                        removeable.append(residue.id)

        # Now that heteroatoms are selected, remove them from the chain
            for residue in removeable:
                chain.detach_child(residue)
            chains.append(chain)

        # Finally, obtain the alpha carbon chain and store it
        # Check if the length of the polypeptide chain is long enough to not be considered a ligand/cofactor (not for nucleotides)
            if len(chain)<=25 and len(next(chain.get_residues()).get_resname())<3:
                structure[0].detach_child(chain.id)
            else:
                alpha_carbon_chains+=1          # Counter for number of structures



        # Check the if we are working with P-Pinteraction or P-Nuc interactions:
        key_chain = [x for x in structure.get_chains()][1]
        chain_type = alpha_carbons_retriever(key_chain, options_verbose)[1]


        if chain_type == "Protein":               # If P-P interaction, we need to have binary interactions (2 CA lists)
            if alpha_carbon_chains!= 2:
                if options_verbose:
                    sys.stderr.write("ERROR: Files do not have right input format to build a complex.\n" )
                continue
            dict_with_PP[id]=structure

        else:
            if chain_type =="DNA":
                if alpha_carbon_chains != 3:       # If DNA-Nuc interaction, we need to have 3 different Seqs in list (Protein, 1st DNA and 2nd DNA).
                    if options_verbose:
                        sys.stderr.write("ERROR: File %s does not have right input format.\n" % (file))
                    continue
                dict_with_dNP[id]=structure

            elif chain_type == "RNA":              # if RNA strands are included, we can have 2 or three chains in the file
                if alpha_carbon_chains == 2:
                    if options_verbose:
                        sys.stdout.write("File %s contains the right input format for a complex including single-strand RNA.\n" % (file))
                    continue
                    dict_with_sNP[id]=structure
                elif alpha_carbon_chains == 3:
                    if options_verbose:
                        sys.stdout.write("File %s contains the right input format for a complex including double-strand RNA.\n" % (file))
                    continue
                    dict_with_dNP[id]=structure

                else:
                    if options_verbose:
                        sys.stderr.write("ERROR: File %s does not have right input format.\n" % (file))
                    continue

    if bool(dict_with_PP) == True:      # Files contain 2x peptide chains
        return (dict_with_PP, "PP")

    elif bool(dict_with_dNP) == True:    # Files contain 2x nucleotide chains and 1x protein chain
        return (dict_with_dNP, "dNP")

    elif bbol(dict_with_sNP) == True:    # Files contain single-strand RNA proteins interaction
        return (dict_with_sNP,"sNP")

    else:
        if options_verbose:
            sys.stdout.write("The provided input files cannot be used to " % (file))

#===================================================================

def alpha_carbons_retriever(chain, options_verbose):
    """
    Get alpha Carbons from input chains (CA for preotein sequence and C4' for DNA/RNA).

    Argument: chain class with the atoms

    Returns: - list of CA or C4 atoms
             - class of molecule we are working with

    """
    nucleic_acids = ['DA','DT','DC','DG','DI','A','U','C','G','I']
    RNA = ['A','U','C','G','I']
    DNA = ['DA','DT','DC','DG','DI']
    atoms = []

    for residue in chain:
        res_type = residue.get_resname().strip()       # Get residue name
        if residue.get_id()[0] == " ":                 # Check if we are dealing with and HET entry

            if res_type not in nucleic_acids:         # If residue is not a nucleic acid
                if 'CA' not in residue:                 # If there are no alpha carbons
                    if options_verbose:                 # And the verboes option has been set in the function: print a message
                                                        # informing about not having CA
                        sys.stdout.write("This residue %d %s doest not have an alpha carbon" % (residue.get_id()[1], res_type))

                else:                                   # If there are alfa cabrons
                    atoms.append(residue['CA'])
                    molecule='Protein'

            elif res_type in DNA:                       # Otherwise, if the residue is a deoxynucleic acid
                molecule = 'DNA'
                atoms.append(residue['C4\''])

            elif res_type in RNA:                       # Finally, if the residue is a nucleic acid
                molecule = 'RNA'
                atoms.append(residue['C4\''])

    return(atoms, molecule)                             # Return the list of alpha carbon atoms and the molecule types

#===================================================================

def dir_path(string):
    """A function to check whether a string is a directory or not"""

    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

#===================================================================

def check_files(path):
    """
    A function to check whether PDB input files have correct format
    """

    work_files=[]
    for file in os.listdir(path):
        if file.endswith(".pdb"):         # Checking if files end with pdb to be consider as input files
            work_files.append(file)


    if not work_files:
        raise ValueError("Check the PDB input files format")
    else:
        os.chdir(path)
        return work_files

#===================================================================

def output_dir(string, options_force):
    """
    A function to check whether outputfile already exists
    """
    if  options_force is False:
        if os.path.isdir(string):
            raise ValueError("ERROR: Directory already exists. Please set -f to True to overwrite the directory")
        else:
            sys.stdout.write("Setting the output directory to %s" % (string))
            os.mkdir(string)

#===================================================================

def align_chains(chain1, chain2):
    """
    Run alignment for two chains of any type
    Return alignment score
    """
    alignment=pairwise2.align.globalxx(chain1,chain2)           # Pairwise alignment of two chains
    alig_score=alignment[0][2]/max(len(chain1),len(chain2))     # calculate alignment score

    return alig_score


#===================================================================

def align_chains_peptides(chain1,chain2):
    """
    A function aligning two chains with each other
    Returns the final alignment score of both peptidic chains
    """
    alpha_carbons=CaPPBuilder()

    chain1_carbons=alpha_carbons.build_peptides(chain1)
    chain1_carbons=chain1_carbons[0].get_sequence()

    chain2_carbons=alpha_carbons.build_peptides(chain2)
    chain2_carbons=chain2_carbons[0].get_sequence()

    return align_chains(chain1_carbons, chain2_carbons)

#===================================================================


def superimpose_chains(ref_structure,alt_structure,threshold, options_verbose):
    """
    Core function to firstly align peptide chains from reference and alternative model.
    Secondly, for those chains found to be similar, superimpose them and obtain
    a dictionary with all the possible superimposition of the chains from the
    two structures (if the superimposition is below a certain RMSD threshold)
    """

    # Initialize variables
    superimpositions={}                                 # dictionary to store possible superimpositions
    best_RMSD=""
    ref_chains=[x for x in ref_structure.get_chains()]  # steady reference chain
    alt_chains=[x for x in alt_structure.get_chains()]  # flexible alternative chain
    sup=Superimposer()   # Superimposer from Biopython

    for ref_chain in ref_chains:
        for alt_chain in alt_chains:
            ref_atoms, ref_molecule = alpha_carbons_retriever(ref_chain, options_verbose)   # get alpha C for ref
            alt_atoms, alt_molecule = alpha_carbons_retriever(alt_chain, options_verbose)   # get alpha C for alt
            if alt_molecule != ref_molecule:
                continue
            if ref_molecule == "Protein":                       # If chain is Protein
                if not align_chains_peptides(ref_chain,alt_chain) > 0.95:   # For the similar chains
                    continue
            elif ref_molecule == "DNA":                         # If chain is DNA append to chain sequence the new ref chain
                ref_chain_seq+=(''.join([x.get_resname()[2] for x in ref_chain]))
                alt_chain_seq+=(''.join([x.get_resname()[2] for x in alt_chain]))
                if not align_chains(ref_chain_seq,alt_chain_seq)>0.95:
                    continue
            elif ref_molecule == "RNA":                         # If chain is RNA append to chain sequence the new ref chain
                ref_chain_seq+=(''.join([x.get_resname()[1] for x in ref_chain]))
                ref_chain_seq+=(''.join([x.get_resname()[1] for x in alt_chain]))
                if not align_chains(ref_chain_seq,alt_chain_seq)>0.95:
                    continue

            sup.set_atoms(ref_atoms,alt_atoms)                  # get rotation and translation matrix
            RMSD=sup.rms                                        # get RMSD for superimposition
            if RMSD < threshold:
                if not best_RMSD or RMSD < best_RMSD:
                    best_RMSD=RMSD
                superimpositions[(ref_chain.id,alt_chain.id)]=sup

    if bool(superimpositions) == True:                                              # If we are able to superimpose any chain
        superimpositions=sorted(superimpositions.items(), key=lambda x:x[1].rms)    # Sort the superimpositions by RMSD
        return (superimpositions,best_RMSD)

#===================================================================

def create_ID(IDs_present):
    """
    Create new ID to make sure it is a non-taken ID
    Input: list of IDs already occupied
    Return: new ID
    """

    Up = list(string.ascii_uppercase)
    Low = list(string.ascii_lowercase)
    Dig = list(string.digits)
    possibilities = set(Up+Low+Dig)              # Set of all acceptable IDs that are possible

    if len(IDs_present)<62:
        possibilities.difference_update(set(IDs_present)) # Update possibilities set by substracting taken ID
        return list(possibilities)[0]

    elif len(IDs_present)>=62:                  # As soon as all possibilities from the set are taken
        for character in possibilities:
            for character2 in possibilities:
                ID = character + character2     # Combine letters to createe new ID
                if ID not in IDs_present:       # Test if new ID already taken
                    return
                else:
                    continue
#========================================================================
def save_structure(structure, options_output, options_verbose, options_force):
    """
    Function saving the generated complexes as pdb files to a folder at the
    output directory provided by the user
    If  a folder already exisst at that directory and the force option is selected,
    it will forcefully overwrite this folder
    """
    if len(list(structure.get_atoms())) > 99999 or len(list(structure.get_chains())) > 62:		# Checks that the structure has has less atoms than the maximum for a PDB, 99,999
        io=MMCIFIO()
        final_complex="final_complex.cif"
    else:
        io=PDBIO()
        final_complex="final_complex.pdb"
    io.set_structure(structure)

    # checking for the existance of directory
    if os.path.exists(options_output):
        if options_force == True:
            shutil.rmtree(options_output)
        else:
            raise ValueError("ERROR: Directory already exists. If you would like to replace it, please provide -f option")
    os.mkdir(options_output)
    os.chdir(options_output)
    # os.mkdir("structure")
    io.save(final_complex)
    if options_verbose:
            sys.stdout.write("YUHU!! The final complex was saved to 'final_complex.pdb' in %s\n" % (options_output))

#========================================================================
def check_for_clashes(ref_structure, added_chain, options_verbose,clash_treshold, success=False):
    """
    Check for clashes between moving structure and refernce structure
    after they have been superimposed.

    Argument: reference structure and moving structure

    Returns: reference structure (added if number of clashes below threshold)

    """

    ref_atoms=[]
    for chain in ref_structure.get_chains():                         # Get all the atom positions in the current reference structure
        ref_atoms.extend(alpha_carbons_retriever(chain,options_verbose)[0])

    moving_atoms = added_chain.get_atoms()

    Neighbor = NeighborSearch(ref_atoms)                             # Using NeighborSearch from Biopython creating an instance Neighbor
    clashes = 0
    for atom in moving_atoms:                                        # Search for possible clashes between the atoms of the chain we want to
        atoms_clashed = Neighbor.search(atom.coord,3 )               # add and the atoms already in the model
        if len(atoms_clashed) > 0:
            clashes+=len(atoms_clashed)

    if clashes < clash_treshold:                                    # If the clashes do not exceed a certain threshold add the chain to the model

        present=[chain.id for chain in ref_structure.get_chains()]
        if added_chain.id in present:
            added_chain.id= create_ID(present)                      # Create random id so it does not clash with the current chain ids in the PDB structure
        ref_structure[0].add(added_chain)
        success = True                                              # Set boolean because chain was added
        if options_verbose:
            sys.stdout.write("The chain %s was added to the model\n" % (added_chain.id))
    else:
        if options_verbose:
            sys.stdout.write("The chain %s was not added to the model: too many clashes.\n" % (added_chain.id))


    return ref_structure, success
