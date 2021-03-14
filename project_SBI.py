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
                    dest="stechiometry",
                    type=dir_path,
                    action="store",
                    help='Directory containing the appropiate stechiometry file for the molecular complex')
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
    # print(structure_data)

# Process Stoichiometry if provided by User
    if options.stechiometry:
        try:
            stech_dict={}
            stech_file=check_stech(options.stechiometry)
            print(stech_file)
        #     with open(stech_file,"r") as file:
        #         for line in file:
        #             a=line.split(":")
        #             if not a[1].isdigit():
        #                 raise ValueError("Not an appropiate stechiometry file")
        #                 stech_dict[a[0]]=a[1]
        except:
            raise ValueError("No stechiometry file found")

        print(stech_dict)
    # except NotADirectoryError as e:
    #     sys.stderr.write("Stechiometry option does not correspond to an existing directory. Please try again.")
    #     exit()
    # if options.stec:
    #     stech_dict={}
    #     script_dir=os.path.dirname(os.path.realpath('__file__'))
    #     rel_path=options.stec
    #     abs_path=os.path.join(script_dir,rel_path)
        #
        # try:
        #     with open(abs_path,"r") as file:
        #         for line in file:
        #             a=line.split(":")
        #             if not a[1].isdigit():
        #                 raise ValueError("Not an appropiate stechiometry file")
        #             stech_dict[a[0]]=a[1]
        # except:
        #     raise ValueError("No stechiometry file found")
        # print(stech_dict)

#Check the heterodimer structures to find the same chain under different chain ids
    if structure_data["heterodimers"]: # Check if we have heterodimers
        same_chains={}
        alpha_carbons=CaPPBuilder()
        heterodimers=structure_data["heterodimers"]
        heterodimer_list=list(heterodimers.keys()) # List holding the files id like "A_B"
        # heterodimer pairwise interaction in list chains1 and chains2
        for j in range(len(heterodimer_list)):
            chains1=heterodimers[heterodimer_list[j]].get_chains()
            for i in range(j+1,len(heterodimer_list)+1):
                chains2=heterodimers[heterodimer_list[i]].get_chains()
                # access the two chains in each heterodimer structure
                for chain in chains1:
                    chain_alpha = alpha_carbons.build_peptides(chain)
                    chain_alpha = chain_alpha[0].get_sequence()
                    for chain2 in chains2:
                        chain_alpha2 = alpha_carbons.build_peptides(chain)
                        chain_alpha2 = chain_alpha2[0].get_sequence()

                        alignment = pairwise2.align.globalxx(chain_alpha, chain_alpha2)
                        if alignment>0.95:
                            if chain2.id in same_chains:
                                continue
                            else:  # save same sequences in dictionary with 2. chainid as key and first chainid as value
                                same_chains[chain2.id]=chain.id
