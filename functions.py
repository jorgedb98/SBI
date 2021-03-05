def read_pdb_files(files):
    """Given a pdb file, read it, remove the heteroatoms and create a dictionary with the chain ids and the structure

    Input:
        PDB File (files argument) with a pairwise interaction

    Output:
    Dictionary with three elements: Chain ids (2) and the structure """

    structure=pdb_parser.get_structure(files)

    chains_ids=[chain.id for chain in structure[0].get_chains()]


def dir_path(string):
    """A function to check whether a string is a directory or not"""
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def check_files(path):
    """A function to check whether inputfiles have correct format"""
    work_files=[]
    my_pattern=re.compile("\w+_\w+_\w+.pdb*")
    for file in os.listdir(path):
        work_files.append(my_pattern.match(file))
    if not work_files:
#    if my_pattern.match(file) == None:
        raise ValueError("Check the input files format")
    else:
        return work_files

def output_dir(string):
    """A function to check whether outputfile already exists"""
    if  options.force is False:
        if os.path.isdir(string):
            raise ValueError("Directory already exists. Please set -f to True to overwrite the directory")
        else:
            sys.stderr.write("Setting the output directory to %s" % (string))
            os.mkdir(string)
