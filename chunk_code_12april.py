
            ref_structure = structure_data[complex_files.pop(0)]
            dna_chain = list(ref_structure.get_chains())[1]                    # takening the 2nd chain as reference DNA for alignment
            #dna_chain=x.get_resname() for x in dna_chain]                     # getting the residue names
            dna_chain_seq=''.join([x.get_resname()[2] for x in dna_chain])     # get residue names for DNA strand

            for complex2 in complex_files:                               # go through list of proteins in complex_files
                moving_structure=structure_data[complex2]
                dna_chain2=list(moving_structure.get_chains())[1]        # get the dna for the moving structure
                dna_chain2_seq=''.join([x.get_resname()[2] for x in dna_chain2])
                if align_chains(dna_chain_seq, dna_chain2_seq) < 0.8:           # if alignment of the two DNA is below threshold
                    dna_chain2=list(moving_structure.get_chains())[2]  # compare with the other DNA strand
                    dna_chain2_seq=''.join([x.get_resname()[2] for x in dna_chain2])
                    if align_chains(dna_chain_seq,dna_chain2_seq) < 0.8:       # if alignment with both DNA strands below threshold
                        print("should continue")
                        continue                                       # no complex build

                sup=Superimposer()                       # alignment is good enough to superimpose complex and complex2
                dna_atoms, molecule=alpha_carbons_retriever(dna_chain, options.verbose)    # transform both DNA chains to list
                dna_atoms2, molecule2=alpha_carbons_retriever(dna_chain2, options.verbose)
                list(dna_chain2.get_atoms())
                print(len(dna_atoms2),len(dna_atoms))

                if len(dna_atoms)>len(dna_atoms2):
                    #retrieve=re.search(dna_chain2_seq,dna_chain_seq).span()
                    print(dna_chain_seq,dna_chain2_seq)
                    continue
                else:
                    #retrieve=re.search(dna_chain_seq,dna_chain2_seq).span()
                    print(dna_chain_seq,dna_chain2_seq)
                    continue
                print(len(dna_atoms))

                sup.set_atoms(dna_atoms2,dna_atoms)    # retrieve rotation and translation matrix
                RMSD=sup.rms                           # get RMSD for superimposition

                if RMSD < 2:                   # if RMSD is below threshold, apply superimposition
                    added_chain= next(moving_structure.get_chains())
                    sup.apply(added_chain.get_atoms())

                    ref_atoms=[]
                    for chain in ref_structure.get_chains():  #Get all the atom positions in the current reference structure
                        ref_atoms.extend(alpha_carbons_retriever(chain,options.verbose)[0])

                    moving_atoms=added_chain.get_atoms()

## could be function check-clashes (we could include the lines above too )
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
# could be function check-clashes


        structure_list=list(structure_data.keys())
        refid=structure_list.pop(0)
        structure_list.append(refid)
        current_stech={refid:1}

        it_count=0

        # SUPERIMPOSE C-alphas of those CHAINS WITH HIGH ALIGNMENT
        # ref_structure = structure_data[refid]           # get first pair as reference structure
        # ref_dna=list(ref_structure.get_chains())[1]     # in the pdb file: 2nd chain as ref DNA
        # ref_dna=[x.get_resname()[2] for x in ref_dna]
        # ref_dna=''.join(ref_dna)
        nc=0
       dictionary to store protein-DNA complexes


        # new_dna = ref_dna.replace(" ", "")
        # data_splited = re.findall('..',new_dna)

        while (nc <= sum(list(stech_file.values()))):      # as long as we need more chains to fullfill the stechiometry
            moveid = structure_list.pop(0)                 # take id for move structure from all structures
            if moveid not in stech_file:                   # if current mov_id not in stechiometry
                moveid = ""                                # do not use it for complex and go to next id
                continue
            if not moveid in current_stech:                # If the count for the current structure id is not initialised, start it
                current_stech[moveid] = 0
            moving_structure = structure_data[moveid]
            move_dna = list(moving_structure.get_chains())[1]      # get second chain of moving_object as DNA
            move_dna = [x.get_resname()[2] for x in move_dna]      # get residues
            alignment = align_chains(ref_dna, ''.join(move_dna))   # align with reference DNA

            if alignment < 0.75:                                     # if alignment below threshold
                move_dna=list(moving_structure.get_chains())[2]      #
                move_dna=[x.get_resname()[2] for x in move_dna]          # If the two seq do not align
                alignment2=align_chains(ref_dna, ''.join(move_dna))  # compare to the reverse DNA strand

                if alignment2 <0.75:                           # if the second DNA chain has an alignment below threhsold too
                    structure_list.append(moveid)              #
                    continue
            ref_dna=''.join(move_dna)

            nc+=1
