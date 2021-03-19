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
    structure_data=read_pdb_files(work_files, options.verbose)
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
                    stech_file[a[0]]=a[1]
        except NotADirectoryError as e:
            sys.stderr.write("Stechiometry option does not correspond to an existing directory. Please try again.\n")
            exit()


#Check the heterodimer structures to find the same chain under different chain ids
    if "heterodimers" in structure_data.keys(): # Check if we have heterodimers
        # print("are you there?")
        same_chains={}
        alpha_carbons=CaPPBuilder()
        heterodimers=structure_data["heterodimers"]
        heterodimer_list=list(heterodimers.keys()) # List holding the files id like "A_B"
        # print(structure_data)

        # hold heterodimer pairwise interaction in list chains1 and chains2
        for j in range(len(heterodimer_list)-1):
            # print("j:",j)
            #print("do you eter here?")
            chains1=heterodimers[heterodimer_list[j]].get_chains()
            # print(list(chains1))
            chains1=heterodimers[heterodimer_list[j]].get_chains()

            ##################################################


            for i in range(j+1,len(heterodimer_list)):
                #chain.id=chain.id+i
                # access the two chains in each heterodimer structure
                chains1=heterodimers[heterodimer_list[j]].get_chains()
                for chain in chains1:
                    #print(heterodimer_list[j])
                    chain_alpha = alpha_carbons.build_peptides(chain)
                    chain_alpha = chain_alpha[0].get_sequence()
                    chains2=heterodimers[heterodimer_list[i]].get_chains()
                    # print(chains2)
                    for chain2 in chains2:

                        # print(chain, chain2)
                        #print(heterodimer_list[i])
                        chain_alpha2 = alpha_carbons.build_peptides(chain2)
                        chain_alpha2 = chain_alpha2[0].get_sequence()

                        #If the chains share the same id, do not compare them since they should be similar
                        if True: #chain.id == chain2.id:
                            alignment = pairwise2.align.globalxx(chain_alpha, chain_alpha2)
                            alig_score=alignment[0][2]/max(len(chain_alpha),len(chain_alpha2))
                            # print(str(chain.id)+ "\t" + chain2.id + "\t" + str(alig_score))
                            if alig_score > 0.95:

                                if not chain2.id in same_chains:
                                    # save same sequences in dictionary with 2. chainid as key and first chainid as value
                                    same_chains[chain2.id]=chain.id

# #Process every scenario differently
#     sup = Superimposer()
#     print(structure_data)
# #Easiest one: All homodimers and no heterodimers
#     if "homodimers" in structure_data.keys() and not "heterodimers" in structure_data.keys():
#         if options.verbose:
#             sys.stderr.write("Input files contain %d homodimers and no heterodimers\n" % len(structure_data["homodimers"]))
#
#         #take one chain as fixed list of atoms
#         #superimpose another binary interaction to the fixed chain using biopython superimposer
#         homodimer_dict=structure_data["homodimers"]
#         ref_model = homodimer_dict[list(homodimer_dict.keys())[0]]
#
#         homodimer_transformed=transform_to_structure(homodimer_dict,"homodimers")
#
#         macrocomplex=Structure.Structure("macrocomplex")
#         macrocomplex.add(Model.Model(0))
#         ref_model = homodimer_transformed[0] # define reference model (first element)
#         macrocomplex.add(list(ref_model.get_chains())[0])
#         macrocomplex.add(list(ref_model.get_chains())[1])
#
#         for alt_model in homodimer_transformed:
#             ref_atoms = []
#             alt_atoms = []
#             for (ref_chain, alt_chain) in zip(ref_model, alt_model):
#                 for ref_res, alt_res in zip(ref_chain, alt_chain):
#                     if ref_res.resname == alt_res.resname and ref_res.id == alt_res.id:
#                         ref_atoms.append(alt_res['CA'])
#                         alt_atoms.append(alt_res['CA'])
#
#             sup.set_atoms(ref_atoms, alt_atoms)
#
#             if ref_model.id == alt_model.id:
#     	        #Check for self/self get zero RMS, zero translation
#     	        #and identity matrix for the rotation.
#                 assert numpy.abs(sup.rms) < 0.0000001
#                 assert numpy.max(numpy.abs(sup.rotran[1])) < 0.000001
#                 assert numpy.max(numpy.abs(sup.rotran[0]) - numpy.identity(3)) < 0.000001
#             else:
#     	        #Update the structure by moving all the atoms in
#     	        #this model (not just the ones used for the alignment)
#                 sup.apply(alt_model.get_atoms())
#                 alt_model.detach_child(list(alt_model.get_chains())[0].id)
#                 sys.stderr.write("RMSD(first model, model %i) = %0.2f" % (alt_model.id, sup.rms))
#
#
#     # print(same_chains) # Checkpoint: same_chains{} holds the ids to chains from different interactions that have over 95% similarity (same ids like "A:A" not inclued)
#
# # All heterodimers:
#     elif "heterodimers" in structure_data.keys() and not "homodimers" in structure_data.keys():
#         print("it enters here")
#         if options.verbose:
#             sys.stderr.write("Input files contain %d heterodimers and no homodimers\n" % len(structure_data["heterodimers"]))
#
#         heterodimer_dict = structure_data["heterodimers"]
#
#         ## Set the maximum number of chains the complex will have and get those files that can be useful
#         if stech_file: # If there is stechiometry file
#             # max_chains = sum(len(list(stech_file.values()))) # The max number = the sum of stechiometry
#             print("1", heterodimer_dict)
#             # get those pairwise interactions that can be used if stechiometry is provided
#             to_remove=[]
#             for structure in heterodimer_dict.values():
#                 chains = [x.id for x in structure.get_chains()]
#                 if not chains[0] in stech_file.keys() and not chains[1] in stech_file.keys():
#                     if not chains[0] in stech_file.keys() and not chains[1] in same_chains.item():
#                         to_remove.append(structure.id) #del interaction file if neither of the chains can be superimposed
#                 elif not chains[0] in stech_file.keys():
#                     if not chains[0] in same_chains.item():
#                         to_remove.append(structure.id)
#                 elif not chains[1] in stech_file.keys():
#                     if not chains[1] in same_chains.items():
#                         to_remove.append(structure.id)
#
#             for id in to_remove:
#                 del heterodimer_dict[id]
#
#             print("2", str(same_chains))
#
#         else: # If not stechiometry provided
#             max_chains = 2*len(heterodimer_dict.keys()) # The max number = twice the number of heterodimer files
#
#
#         heterodimer_transformed = transform_to_structure(heterodimer_dict, "heterodimer")
#
#         # for heterodimers in len(heterodimer_dict.items()): # iterate over list of structure using each as ref model
#         #     ref_model = heterodimer_transformed[heterodimer]
#         #     i = 0
#         #     while i < max_chains: # Iterate while number of max chains is not achieved
#
#
#  #
# # #Homodimers and heterodimer
# #     else:
# #         if options.verbose:
# #             sys.stderr.write("Input files contain %d homodimers and %d heterodimers\n" % (len(structure_data["homodimers"],len(structure_data["heterodimers"])))


#### TESTING AITORS NEW CODE
    het_list = list(structure_data["heterodimers"].keys())
    heterodimers = structure_data["heterodimers"]

# SUPERIMPOSE CHAINS WITH HIGH ALIGNMENT
    ref_structure = heterodimers[het_list[0]] # get first pair as reference structure
    for i in range(1,len(het_list)):
        moving_structure = heterodimers[het_list[i]]
        superimposition = superimpose_chains(ref_structure, moving_structure,2)

        if bool(superimposition) == False:
            to_end=het_list.pop(0) # pop the first heterodimer pair if no superimposition could be made
            het_list.append(to_end) # Append it to the end of the list to see if it can be superimposed latter
            continue
    #Program continues if there are superimpositions

        for possibility, sup in superimposition[0]: #Iterate over the dictionary with superimpositions, not the RMSD
            added_chain = [chain for chain in moving_structure if chain.id != possibility[1]][0]
            sup.apply(added_chain.get_atoms()) # apply rotation matrix to moving structure
            ref_atoms = list(ref_structure.get_atoms()) # atoms of reference model
            moving_atoms = list(added_chain.get_atoms()) # atoms of moving model
            print(moving_atoms)
            # Neighbor = NeighborSearch(ref_atoms) # using NeighborSearch from Biopython creating an instance Neighbor
            # clashes = []
            # for atoms in moving_atoms:
            #     atoms_clashed = Neighbor.search(atoms.coord,5) # check for atoms overlapping between moving structure and object Neighbor from reference structure
            #     # print(atoms_clashed)
            #     if len(atoms_clashed) > 0:				#if there are clashes
            #         clashes.extend(atoms_clashed)       # add clashing atoms to clashlist
            #
            #
            #     # print(clashes)
            # print(len(clashes))
