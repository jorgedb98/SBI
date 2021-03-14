# Work pipeline
___First draft of a pipeline___

## Step 0:
Read pdb files from the folder

Given:
A2B2  (stechiometry by user)

Files:
file1: A-B
file2: A-C
file3: A2-D2 (where none of the 2 chains are the same as in the previous files)


## Step 1:
pdb parser to get chains for alignment:

Seperate the binary interactions into Homo and Heterodimer (95% similarity)
--> Get Homodimer-dictionary and Heterodimer-dictionary with key:ID and value:Structure




Compose a sequence alignment between all our files
file n to file n+1

-> exclude files, if it does not pass the threshold with any of the other files

Outcome: dictionary holding all files important for our model building
{key: A-B, value: PDB_Structure}

## Step 2:
Superimpose
to build all possible models

DONE :)


# Pipeline of the "Massonix" github
Step 0:
    LIST with all pdb files

    for elements in filelist:
      GetStructures function to get structures
      DICTIONNARY with: filepath(structure)

    __stichiometry__
      they assume stechiometry given by the user


Step 1:
Remove chains if they are only found once (cannot be superimposed)

Step 2: SUPERIMPOSITION
Recursive Function:
  - superimpose two complexes
  if no clashes (backbone overlap) -> take this as new model
    else: stop

    repeat with newly build complex as target

  --> to preserve from unnecessary computational cost:

      a. STOP if complex was build already:
          superimpose the newly build complex to previous once
          if RMSD < 2: this complex was build before

      b. STOP out if stechiometry is not met


Step3:
They have an optional argument for optimization
if given, they run the function "Optimizemodel" on the output_dir_files
Optimization of the model(s) by MODELER (rearranging loops and sidechains to minimize energy)



                    A & B --> 0.95    dictionary key[B]=A
                    B & F --> 0.95    dictionary key[F]=B
                    # A & D --> 0.95    dictionary key[D]=A
                    D & A --> 0.95      dictionary would be key[A]=D
                    D & B
