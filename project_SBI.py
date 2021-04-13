import os, sys, gzip, argparse, re, glob, numpy

from functions import *
from Bio.PDB import *


parser=argparse.ArgumentParser(description="SBI_PYT program\n")
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
                    help="Directory for the output files. If the directory already exists, check the -f parameter.\n")

parser.add_argument('-f','--force',
                    default=False,
                    dest="force",
                    action="store_true",
                    help="In case the output directory exists, set the value to true in order to replace it. Otherwise, the output will be stored in a subdirectory.\n")

parser.add_argument('-v','--verbose',
                    default=False,
                    dest="verbose",
                    action="store_true",
                    help="Set the value to true if you want to console output of the progress\n")

parser.add_argument('-m','--max_iter',
                    default=100,
                    dest="max_iterations",
                    action="store",
                    type=int,
                    help="The number maximum iterations the program will try to expand the structure. \
                        WARNING: If a high number of files is provided, the default value is set to 100, thus it may leave the structure unfinished.\n")

parser.add_argument('-n','--nuc',
                    default=None,
                    dest="nuc",
                    action="store",
                    help="Only needed when nucleotide-protein interactions are provided. The path to the file holding the model nucelic acid pdb file (for DNA or RNA) will be stored.\n")

parser.add_argument('-c','--clashes',
                    default=30,
                    dest="clashes",
                    action='store',
                    type=int,
                    help='Maximum number of clashes allowed when trying to superimpose two structures.\n')

parser.add_argument('-t','--threshold',
                    default=3,
                    type=int,
                    action='store',
                    dest='threshold',
                    help='Max RMSD treshold allowed when superimposing two structures.\n')

parser.add_ar

options=parser.parse_args()

if __name__=="__main__":
    absolute_path = os.getcwd()
    try:
        work_files=check_files(options.input)
    except NotADirectoryError as e:
        sys.stderr.write("Input option does not correspond to an existing directory. Please try again.\n")
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
            sys.stderr.write("You have not provided a stoichiometry. Your model will be built using the default value\n")
        stech_file={}
        if interaction =="PP":
            for key in structure_data:
                stech_file[key]=1
        else:
            keys=set([x.split('.')[0] for x in structure_data.keys()])
            for a in keys:
                stech_file[a]=1
    interaction ='PP'
    if interaction == "PP": # When files contain PP complex
        if options.verbose:
            sys.stderr.write("The files provided contain a Protein-Protein interaction.\n")

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
                if options.verbose:
                    sys.stderr.write("No superimpositions found structure %s.\n" %(moveid) )

                it_count+=1                                   # increase count of iterations if no iteration was found
                prot_list.append(moveid)                      # if no superimposition found, bring current moving structure to end of list
                continue

            # Program continues if there are superimpositions
            for possibility, sup in superimposition[0]: #Iterate over the dictionary with superimpositions, not the RMSD
                added_chain = [chain for chain in moving_structure.get_chains() if chain.id != possibility[1]][0]
                sup.apply(added_chain.get_atoms()) # apply rotation matrix to moving structure

                ref_structure, success =check_for_clashes(ref_structure, added_chain, options.verbose, options.clashes)

                if success is True:
                    nc+=1

            if current_stech[moveid] != stech_file[moveid]:   # If structure not as in stechiometry
                prot_list.append(moveid)                      # Append it to the end of the list to see if it can be superimposed later

            it_count+=1
            if it_count == options.max_iterations:
                break
                    #Save the structure
                    # If number of ids taken is lower or eq to 62

        save_structure(ref_structure[0], options.output, options.verbose, options.force)


    elif interaction == "dNP":          # when file contains nucleotide chains and protein chain
        if options.verbose:
            sys.stderr.write("The files provided contain a double-strand Nucleotide-Protein interaction.\n")

        if not options.nuc:             # User did not provide the nucleotide sequence needed for building a complex
            raise ValueError("It seems you are trying to build a macrocomplex for nucleic acid-protein interacion but you are lacking the reference nucleic acid structure. Please, provide it.\n")
        id="ref_dna"

        pdb_parser=PDBParser(PERMISSIVE=1, QUIET=True)

        ref_dna = pdb_parser.get_structure(id,options.nuc)    # read nucelotide strand from user input as reference
        ref_dna_chains = list(ref_dna.get_chains())           # get the strands of reference nucleotide chain

        ref_atoms,molecule=alpha_carbons_retriever(ref_dna_chains[0], options.verbose)
        ref_atoms2, molecule=alpha_carbons_retriever(ref_dna_chains[1], options.verbose)
        ref_atoms.extend(ref_atoms2)

        if molecule =="DNA":
            ref_chain_seq = ''.join([x.get_resname()[2] for x in ref_dna_chains[0]]) # get sequence in case of DNA
            ref_chain_seq+=(''.join([x.get_resname()[2] for x in ref_dna_chains[1]]))
        else:
            ref_chain_seq = ''.join([x.get_resname()[1] for x in ref_dna_chains[0]]) # get sequence in case of RNA

        # print(ref_chain_seq)
        for complex in stech_file.keys():                     # for all complexes that will be need for the stechiometry
            if options.verbose:
                sys.stderr.write("\nProcessing complex %s\n" % (complex))
            complex_files = [x for x in structure_data.keys() if x.split('.')[0] == complex]   # if the first name of the complex_id is equal store it to list complex_files
            # i=0
            # while(i<stech_file[complex]):               # while the complex is needed as many times as indicated in the stechiometry
            a=0
            for protein_interaction in complex_files:
                i=0
                # print(protein_interaction)
                # print(list(structure_data[protein_interaction].get_chains()))
                dna_chain1 = list(structure_data[protein_interaction].get_chains())[1]      # get first nucleotide chain in structure
                dna_chain_seq = ''.join([x.get_resname()[2] for x in dna_chain1])           # get sequence
                #print(re.match(dna_chain_seq,ref_chain_seq))
                possible_locations=[]
                possible_locations=[x.span() for x in re.finditer(dna_chain_seq, ref_chain_seq)]    # span of matches
                #print(possible_locations)

                if not possible_locations:
                    dna_chain1 = list(structure_data[protein_interaction].get_chains())[2]      # get first nucleotide chain in structure
                    dna_chain_seq = ''.join([x.get_resname()[2] for x in dna_chain1])           # get sequence
                    #print("Looking for a pattern as longest as %d" %(len(dna_chain_seq)))
                    possible_locations=[x.span() for x in re.finditer(dna_chain_seq, ref_chain_seq)]

                    if not possible_locations:
                        if options.verbose:
                            sys.stderr.write("For chain %s no pattern matched with the reference nucleotide sequence\n" % (dna_chain1.id))
                        continue

                while(i < stech_file[complex]):                     # while the complex is needed as many times as indicated in the stechiometry
                    start=possible_locations[i][0]          #Set the coordinate of where the first base of our DNA fragment is located
                    #print(start)
                    end=possible_locations[i][1]            #Set the coordinate of where the last base of our DNA fragment is located
                    #print(end)
                    possible_locations.append(possible_locations[i])  #Append the coordinates to the end of the list, in case we may need to reuse them
                    tmp_ref_atoms=ref_atoms[start:end]      #Extract the atoms our atoms from the big reference DNA from one of the possible locations(if available)

                    moving_atoms, molecule=alpha_carbons_retriever(dna_chain1, options.verbose)

                    #print(moving_atoms,tmp_ref_atoms)
                    if len(moving_atoms) != len(tmp_ref_atoms):
                        if options.verbose:
                            sys.stderr.write("Lengths are different. Chain %s will be ignored.\n" % (dna_chain1.id))
                        i += 1
                        continue
                        #print(moving_atoms)
                        #print(list(dna_chain1.get_atoms()))
                    sup=Superimposer()
                    sup.set_atoms(tmp_ref_atoms, moving_atoms)
                    # print(dna_chain_seq, ref_chain_seq[start:end])

                    if not sup.rms < 3:     # Check if RMSD is below threshold
                        i += 1
                        if options.verbose:
                            sys.stderr.write("RMSD score was %d. Chain %s will be ignored.\n" % (sup.rms, dna_chain1.id))
                        continue

                    # Apply rotation matrix to protein chain of structure:
                    chain_to_add=list(structure_data[protein_interaction].get_chains())[0]  # list of proteins for superimposed DNA strand
                    sup.apply(chain_to_add.get_atoms())     # apply translation and rotation matrix to chain

                    ref_dna,success=check_for_clashes(ref_structure=ref_dna, added_chain=chain_to_add, clash_treshold=options.clashes, options_verbose=options.verbose)
                    if success is True:
                        if options.verbose:
                            sys.stderr.write("RMSD score was %d. Chain %s was added.\n" % (sup.rms, dna_chain1.id))
                    #     i+=1 # increase stoichometry counter since new chain was added
                    i=+1

        save_structure(ref_dna[0], options.output, options.verbose, options.force)
        #else cannot save as pdb -> save as MMCIFIO

    elif interaction == "sNP":
        if options.verbose:
            sys.stderr.write("The files provided contain a single-strand RNA  Protein interaction.\n")

        if not options.nuc:             # User did not provide the nucleotide sequence needed for building a complex
            raise ValueError("It seems you are trying to build a macrocomplex for nucleic acid-protein interacion but you are lacking the reference nucleic acid structure. Please, provide it.\n")


    else:
        sys.stderr.write("We are so sorry to tell you your files don't have Protein-Protein nor Nucleotide-Protein interactions :(")
