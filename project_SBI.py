import os, sys, gzip, argparse, re, glob, numpy

from functions import *
from Bio.PDB import *

parser=argparse.ArgumentParser(description="SBI_PYT program")
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-i','--input',
                    required=True,
                    dest="input",
                    action="store",
                    help='file directory. The files that will be sought for must have the format .pdb. ',
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

parser.add_argument('-n','--nuc',
                    default=None,
                    dest="nuc",
                    action="store",
                    help="Only needed when nucleotide-protein interactions are provided. The path to the file holding the model nucelic acid pdb file (for DNA or RNA) will be stored.")

options=parser.parse_args()


if __name__=="__main__":
    absolute_path = os.getcwd()
    try:
        work_files=check_files(options.input)
    except NotADirectoryError as e:
        sys.stderr.write("Input option does not correspond to an existing directory. Please try again.")
        exit()

    if options.verbose:
        sys.stderr.write("%d files found in %s \n" % (len(work_files),options.input))

# Get all the files as structures in a dictionary with the pdb parser from the Bio package and the function we built around it
    structure_data, interaction=read_pdb_files(work_files, options.verbose)
    os.chdir(absolute_path) # After processing input files get back to main working directory

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
    #print(stech_file)

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
        while(nc<=sum(list(stech_file.values()))):    # Iterate while number of chains is lower than SUM of number of chains in stech file
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

    # could be function check-clashes (we could include the lines above too )
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
        # could be function check-clashes end

            if current_stech[moveid] != stech_file[moveid]:   # If structure not as in stechiometry
                prot_list.append(moveid)                      # Append it to the end of the list to see if it can be superimposed later

            it_count+=1
            if it_count == options.max_iterations:
                break
                    #Save the structure
                    # If number of ids taken is lower or eq to 62

        save_structure(ref_structure[0], options.output, options.verbose, options.force)


    elif interaction == "dNP":          # when file contains DNA chains and protein chain
        if options.verbose:
            sys.stderr.write("The files provided contain a double-strand Nucleotide-Protein interaction.")

        if not options.nuc:             # User did not provide the nucleotide sequence needed for building a complex
            raise ValueError("It seems you are trying to build a macrocomplex for nucleic acid-protein interacion but you are lacking the reference nucleic acid structure. Please, provide it.")
        id="ref_dna"

        pdb_parser=PDBParser(PERMISSIVE=1, QUIET=True)

        ref_dna = pdb_parser.get_structure(id,options.nuc)    # read nucelotide strand from user input as reference
        ref_dna_chains = list(ref_dna.get_chains())           # get the strands of reference nucleotide chain
        ref_atoms,molecule=alpha_carbons_retriever(ref_dna_chains[0], options.verbose)
        if molecule =="DNA":
            ref_chain_seq = ''.join([x.get_resname()[1] for x in ref_dna_chains[0]]) # get sequence in case of DNA
        else:
            ref_chain_seq = ''.join([x.get_resname()[0] for x in ref_dna_chains[0]]) # get sequence in case of RNA

        print(ref_chain_seq)
        macrocomplex_dict = {}
        for complex in stech_file.keys():                     # for all complexes that will be need for the stechiometry
            complex_id=complex.split('.')[0]
            complex_files = [x for x in structure_data.keys() if x.split('.')[0] == complex_id]   # if the first name of the complex_id is equal store it to list complex_files
            i=0
            while(i<stech_file[complex]):                     # while the complex is needed as many times as indicated in the stechiometry
                for protein_interaction in complex_files:
                    print(list(structure_data[protein_interaction].get_chains()))
                    dna_chain1 = list(structure_data[protein_interaction].get_chains())[1]      # get first nucleotide chain in structure
                    dna_chain_seq = ''.join([x.get_resname()[2] for x in dna_chain1])           # get sequence
                    # print(dna_chain_seq)
                    print(re.match(dna_chain_seq,ref_chain_seq))
                    possible_locations=[x.span() for x in re.finditer(dna_chain_seq, ref_chain_seq)]    # span of matches
                    #print(len(possible_locations))
                    #print(dna_chain_seq, ref_chain_seq)
                i += 1

            # ref_structure = structure_data[complex_files.pop(0)]
            # dna_chain = list(ref_structure.get_chains())[1]                    # takening the 2nd chain as reference DNA for alignment
            # #dna_chain=x.get_resname() for x in dna_chain]                     # getting the residue names
            # dna_chain_seq=''.join([x.get_resname()[2] for x in dna_chain])     # get residue names for DNA strand

            # for complex2 in complex_files:                               # go through list of proteins in complex_files
            #     moving_structure=structure_data[complex2]
            #     dna_chain2=list(moving_structure.get_chains())[1]        # get the dna for the moving structure
            #     dna_chain2_seq=''.join([x.get_resname()[2] for x in dna_chain2])
            #     if align_chains(dna_chain_seq, dna_chain2_seq) < 0.8:           # if alignment of the two DNA is below threshold
            #         dna_chain2=list(moving_structure.get_chains())[2]  # compare with the other DNA strand
            #         dna_chain2_seq=''.join([x.get_resname()[2] for x in dna_chain2])
            #         if align_chains(dna_chain_seq,dna_chain2_seq) < 0.8:       # if alignment with both DNA strands below threshold
            #             print("should continue")
            #             continue                                       # no complex build
            #
            #     sup=Superimposer()                       # alignment is good enough to superimpose complex and complex2
            #     dna_atoms, molecule=alpha_carbons_retriever(dna_chain, options.verbose)    # transform both DNA chains to list
            #     dna_atoms2, molecule2=alpha_carbons_retriever(dna_chain2, options.verbose)
                # list(dna_chain2.get_atoms())
                # print(len(dna_atoms2),len(dna_atoms))

                # if len(dna_atoms)>len(dna_atoms2):
                #     #retrieve=re.search(dna_chain2_seq,dna_chain_seq).span()
                #     print(dna_chain_seq,dna_chain2_seq)
                #     continue
                # else:
                #     #retrieve=re.search(dna_chain_seq,dna_chain2_seq).span()
                #     print(dna_chain_seq,dna_chain2_seq)
                #     continue
                # print(len(dna_atoms))

#                 sup.set_atoms(dna_atoms2,dna_atoms)    # retrieve rotation and translation matrix
#                 RMSD=sup.rms                           # get RMSD for superimposition
#
#                 if RMSD < 2:                   # if RMSD is below threshold, apply superimposition
#                     added_chain= next(moving_structure.get_chains())
#                     sup.apply(added_chain.get_atoms())
#
#                     ref_atoms=[]
#                     for chain in ref_structure.get_chains():  #Get all the atom positions in the current reference structure
#                         ref_atoms.extend(alpha_carbons_retriever(chain,options.verbose)[0])
#
#                     moving_atoms=added_chain.get_atoms()
#
# ## could be function check-clashes (we could include the lines above too )
#                     Neighbor = NeighborSearch(ref_atoms) # using NeighborSearch from Biopython creating an instance Neighbor
#                     clashes = 0
#                     for atom in moving_atoms: #Search for possible clashes between the atoms of the chain we want to add and the atoms already in the model
#                         atoms_clashed = Neighbor.search(atom.coord,5)
#
#                         if len(atoms_clashed) > 0:
#                             clashes+=len(atoms_clashed)
#
#                     if clashes < 30:   #If the clashes do not exceed a certain threshold add the chain to the model
#                         present=[chain.id for chain in ref_structure.get_chains()]
#                         if added_chain.id in present:
#                             added_chain.id= create_ID(present) #Change the id so it does not clash with the current chain ids in the PDB structure
#                         ref_structure[0].add(added_chain)
## could be function check-clashes


        #structure_list=list(structure_data.keys())
        #refid=structure_list.pop(0)
        #structure_list.append(refid)
        #current_stech={refid:1}

        # it_count=0
        #
        # # SUPERIMPOSE C-alphas of those CHAINS WITH HIGH ALIGNMENT
        # # ref_structure = structure_data[refid]           # get first pair as reference structure
        # # ref_dna=list(ref_structure.get_chains())[1]     # in the pdb file: 2nd chain as ref DNA
        # # ref_dna=[x.get_resname()[2] for x in ref_dna]
        # # ref_dna=''.join(ref_dna)
        # nc=0
       # dictionary to store protein-DNA complexes


        # # new_dna = ref_dna.replace(" ", "")
        # # data_splited = re.findall('..',new_dna)
        #
        # while (nc <= sum(list(stech_file.values()))):      # as long as we need more chains to fullfill the stechiometry
        #     moveid = structure_list.pop(0)                 # take id for move structure from all structures
        #     if moveid not in stech_file:                   # if current mov_id not in stechiometry
        #         moveid = ""                                # do not use it for complex and go to next id
        #         continue
        #     if not moveid in current_stech:                # If the count for the current structure id is not initialised, start it
        #         current_stech[moveid] = 0
        #     moving_structure = structure_data[moveid]
        #     move_dna = list(moving_structure.get_chains())[1]      # get second chain of moving_object as DNA
        #     move_dna = [x.get_resname()[2] for x in move_dna]      # get residues
        #     alignment = align_chains(ref_dna, ''.join(move_dna))   # align with reference DNA
        #
        #     if alignment < 0.75:                                     # if alignment below threshold
        #         move_dna=list(moving_structure.get_chains())[2]      #
        #         move_dna=[x.get_resname()[2] for x in move_dna]          # If the two seq do not align
        #         alignment2=align_chains(ref_dna, ''.join(move_dna))  # compare to the reverse DNA strand
        #
        #         if alignment2 <0.75:                           # if the second DNA chain has an alignment below threhsold too
        #             structure_list.append(moveid)              #
        #             continue
        #     ref_dna=''.join(move_dna)
        #
        #     nc+=1

        save_structure(ref_structure[0], options.output, options.verbose, options.force)
        #else cannot save as pdb -> save as MMCIFIO

    elif interaction == "sNP":
        if options.verbose:
            sys.stderr.write("The files provided contain a single-strand RNA  Protein interaction.")

        if not options.nuc:             # User did not provide the nucleotide sequence needed for building a complex
            raise ValueError("It seems you are trying to build a macrocomplex for nucleic acid-protein interacion but you are lacking the reference nucleic acid structure. Please, provide it.")


    else:
        sys.stderr.write("We are so sorry to tell you your files don't have Protein-Protein nor Nucleotide-Protein interactions :(")
