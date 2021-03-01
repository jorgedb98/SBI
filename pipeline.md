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
pdb parser to get chains for alignment

Compose a sequence alignment between all our files
file n to file n+1

-> exclude files, if it does not pass the threshold with any of the other files

Outcome: dictionary holding all files important for our model building
{key: A-B, value: PDB_Structure}

## Step 2:
Superimpose
to build all possible models

DONE :)
