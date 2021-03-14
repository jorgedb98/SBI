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
                    default=None,
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
        sys.stderr.write("%d files found in $s \n" % (len(work_files),options.input))

#Get all the files as structures in a dictionary with the pdb parser from the Bio package and the function we built around it
    structure_data=read_pdb_files(work_files)
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
        same_chains={}
        alpha_carbons=CaPPBuilder()
        heterodimers=structure_data["heterodimers"]
        heterodimer_list=list(heterodimers.keys()) # List holding the files id like "A_B"

        # hold heterodimer pairwise interaction in list chains1 and chains2
        for j in range(len(heterodimer_list)-1):
            chains1=heterodimers[heterodimer_list[j]].get_chains()

            for i in range(j+1,len(heterodimer_list)+1):
                # access the two chains in each heterodimer structure
                for chain in chains1:
                    chain_alpha = alpha_carbons.build_peptides(chain)
                    chain_alpha = chain_alpha[0].get_sequence()
                    chains2=heterodimers[heterodimer_list[i]].get_chains()

                    for chain2 in chains2:
                        chain_alpha2 = alpha_carbons.build_peptides(chain2)
                        chain_alpha2 = chain_alpha2[0].get_sequence()

                        #If the chains share the same id, do not compare them since they should be similar
                        if chain.id == chain2.id:
                            continue

                        alignment = pairwise2.align.globalxx(chain_alpha, chain_alpha2)
                        alig_score=alignment[0][2]/max(len(chain_alpha),len(chain_alpha2))
                        if alig_score > 0.95:
                            if chain2.id in same_chains:
                                pass
                            else:  # save same sequences in dictionary with 2. chainid as key and first chainid as value
                                same_chains[chain2.id]=chain.id
    else:
        pass

#Process every scenario differently
#Easiest one: All homodimers and no heterodimers
    if "homodimers" in structure_data.keys() and not "heterodimers" in structure_data.keys():
        if options.verbose:
            sys.stderr.write("Input files contain %d homodimers and no heterodimers\n" % len(structure_data["homodimers"]))

        #take one chain as fixed list of atoms
        #superimpose another binary interaction to the fixed chain using biopython superimposer
        homodimer_dict=structure_data["homodimers"]
        ref_model = homodimer_dict[list(homodimer_dict.keys())[0]]

        sup = Superimposer()

        homodimer_transformed=transform_to_structure(homodimer_dict,"homodimers")

        macrocomplex=Structure.Structure("macrocomplex")
        macrocomplex.add(Model.Model(0))
        ref_model=homodimer_transformed[0] # define reference model (first element)
        macrocomplex.add(list(ref_model.get_chains())[0])
        macrocomplex.add(list(ref_model.get_chains())[1])

        for alt_model in homodimer_transformed:
            ref_atoms = []
            alt_atoms = []
            for (ref_chain, alt_chain) in zip(ref_model, alt_model):
                for ref_res, alt_res in zip(ref_chain, alt_chain):
                    if ref_res.resname == alt_res.resname and ref_res.id == alt_res.id:
                        ref_atoms.append(alt_res['CA'])
                        alt_atoms.append(alt_res['CA'])

            sup.set_atoms(ref_atoms, alt_atoms)

            if ref_model.id == alt_model.id:
    	        #Check for self/self get zero RMS, zero translation
    	        #and identity matrix for the rotation.
                assert numpy.abs(sup.rms) < 0.0000001
                assert numpy.max(numpy.abs(sup.rotran[1])) < 0.000001
                assert numpy.max(numpy.abs(sup.rotran[0]) - numpy.identity(3)) < 0.000001
            else:
    	        #Update the structure by moving all the atoms in
    	        #this model (not just the ones used for the alignment)
                sup.apply(alt_model.get_atoms())
                alt_model.detach_child(list(alt_model.get_chains())[0].id)
                sys.stderr.write("RMS(first model, model %i) = %0.2f" % (alt_model.id, sup.rms))



    # print(same_chains) # Checkpoint: same_chains{} holds the ids to chains from different interactions that have over 95% similarity (same ids like "A:A" not inclued)

#All heterodimers:
# elif bool(structure_data["heterodimers"]) and not bool(structure_data["homodimers"]):
#         if options.verbose:
#             sys.stderr.write("Input files contain %d heterodimers and no homodimers\n" % len(structure_data["heterodimers"]))
#
# #Homodimers and heterodimers
#     else:
#         if options.verbose:
#             sys.stderr.write("Input files contain %d homodimers and %d heterodimers\n" % (len(structure_data["homodimers"],len(structure_data["heterodimers"])))
