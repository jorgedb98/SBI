import os, sys, gzip, argparse, re, glob, numpy

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
                    dest="stechiometry",
                    action="store",
                    help='Path to the file where molecular stechiometry for the final model is stored.')
                    #stoichiometry file provided by user in fomrat: identifier ":" counter

requiredNamed.add_argument('-o','--output',
                    required=True,
                    dest="output",
                    action="store",
                    type=str,
                    help="Directory for the output files. If the directory already exists, check the -f parameter.")

parser.add_argument('-f','--force',
                    default=False,
                    dest="force",
                    action="store_true",
                    help="In case the output directory exists, set the value to true in order to replace it. Otherwise, the output will be stored in a subdirectory.")

parser.add_argument('-v','--verbose',
                    default=False,
                    dest="verbose",
                    action="store_true",
                    help="Set the value to true if you want to console output of the progress")

parser.add_argument('-m','--max_iter',
                    default=100,
                    dest="max_iterations",
                    action="store",
                    type=int,
                    help="The number maximum iterations the program will try to expand the structure. \
                        WARNING: If a high number of files is provided, the default value is set to 100, thus it may leave the structure unfinished.")
options=parser.parse_args()


if __name__=="__main__":

    try:
        work_files=check_files(options.input)
    except NotADirectoryError as e:
        sys.stderr.write("Input option does not correspond to an existing directory. Please try again.")
        exit()

    if options.verbose:
        sys.stderr.write("%d files found in %s \n" % (len(work_files),options.input))

#Get all the files as structures in a dictionary with the pdb parser from the Bio package and the function we built around it
    structure_data, interaction=read_pdb_files(work_files, options.verbose)
    os.chdir("../") # after processing input files get back to main working directory

# Process Stoichiometry if provided by User
    if options.stechiometry:
        if options.verbose:
            sys.stderr.write("You have provided the stoichiometry found in %s\n" % (options.stechiometry))
        try:
            stech_file={}
            with open(options.stechiometry, "r") as file:
                for line in file:
                    a=line.strip().split(":")
                    stech_file[a[0]]=int(a[1])
        except NotADirectoryError as e:
            sys.stderr.write("Stechiometry option does not correspond to an existing directory. Please try again.\n")
            exit()
    else:
        if options.verbose:
            sys.stderr.write("You have not provided a stoichiometry. Your model will be built using the default value")
        stech_file={}
        for key in structure_data:
            stech_file[key]=1
    print(stech_file)

    if interaction == "PP": # When files contain PP complex
        if options.verbose:
            sys.stderr.write("The files provided contain a Protein-Protein interaction.")

        prot_list = list(structure_data.keys())
        refid=prot_list.pop(0)
        prot_list.append(refid)
        current_stech={refid:1}

        it_count=0

        # SUPERIMPOSE C-alphas of those CHAINS WITH HIGH ALIGNMENT
        ref_structure = structure_data[refid] # Get first pair as reference structure
        nc=2
        while(nc<=sum(list(stech_file.values()))):    # Iterate while number of chains is lower than sum of number of chains in stech file
            moveid=prot_list.pop(0)
            if moveid not in stech_file:
                moveid=""
                continue
            if not moveid in current_stech:                 # If the count for the current structure id is not initialised, start it
                current_stech[moveid]=0
            moving_structure = structure_data[moveid]
            superimposition = superimpose_chains(ref_structure, moving_structure, 2, options.verbose)
            current_stech[moveid]+=1
            nc +=1
            if current_stech[moveid] != stech_file[moveid]:   # If structure not as in stechiometry
                prot_list.append(moveid)                      # Append it to the end of the list to see if it can be superimposed later
            if bool(superimposition) == False:                # If no superimposition was made, continue to next structure
                it_count+=1                                   # increase count of iterations if no iteration was found
                prot_list.append(moveid)                      # if no superimposition found, bring current moving structure to end of list
                continue

        #Program continues if there are superimpositions
            for possibility, sup in superimposition[0]: #Iterate over the dictionary with superimpositions, not the RMSD
                added_chain = [chain for chain in moving_structure.get_chains() if chain.id != possibility[1]][0]
                sup.apply(added_chain.get_atoms()) # apply rotation matrix to moving structure

                moving_atoms, moving_molecule = alpha_carbons_retriever(added_chain, options.verbose)

                ref_atoms=[]
                for chain in ref_structure.get_chains():  #Get all the atom positions in the current reference structure
                    ref_atoms.extend(alpha_carbons_retriever(chain,options.verbose)[0])

                Neighbor = NeighborSearch(ref_atoms) # using NeighborSearch from Biopython creating an instance Neighbor
                clashes = 0
                for atom in moving_atoms: #Search for possible clashes between the atoms of the chain we want to add and the atoms already in the model
                    atoms_clashed = Neighbor.search(atom.coord,5)

                    if len(atoms_clashed) > 0:
                        clashes+=len(atoms_clashed)

                if clashes < 30:   #If the clashes do not exceed a certain threshold add the chain to the model
                    present=[chain.id for chain in ref_structure.get_chains()]
                    if added_chain.id in present:
                        added_chain.id= create_ID(present) #Change the id so it does not clash with the current chain ids in the PDB structure
                    ref_structure[0].add(added_chain)
                    nc+=1

                it_count+=1
                if it_count == options.max_iterations:
                    break
                    #Save the structure
                    # If number of ids taken is lower or eq to 62
        io=PDBIO()
        io.set_structure(ref_structure[0])
        io.save("test_ourstc.pdb")
        #else cannot save as pdb -> save as MMCIFIO

    elif interaction == "NP":  # when file contains NP complex
        if options.verbose:
            sys.stderr.write("The files provided contain a Nucleotide-Protein interaction.")


    else:
        sys.stderr.write("We are so sorry to tell you your files don't have Protein-Protein nor Nucleotide-Protein interactions :(")
