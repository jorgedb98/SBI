import os, sys, gzip, argparse, re, glob

def dir_path(string):
    """A function to check whether a string is a directory or not"""
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def check_files(path):
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
    if  options.force is False:
        if os.path.isdir(string):
            raise ValueError("Directory already exists. Please set -f to True to overwrite the directory")
        else:
            sys.stderr.write("Setting the output directory to %s" % (string))
            os.mkdir(string)

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
                    dest="stec",
                    type=dir_path,
                    action="store",
                    help='Directory containing the appropiate stechiometry file for the molecular complex')

requiredNamed.add_argument('-o','--output',
                    dest="output",
                    required=True,
                    action="store",
                    type=str,
                    help="Directory for the output files. If the directory already exists, check the -f parameter.")

parser.add_argument('-f','--force',
                    dest="force",
                    default=False,
                    action="store_true",
                    help="In case the output directory exists, set the value to true in order to replace it. Otherwise, the output will be stored in a subdirectory.")

parser.add_argument('-v','--verbose',
                    dest="verbose",
                    default=None,
                    action="store_true",
                    help="Set the value to true if you want to console output of the progress")

options=parser.parse_args()

if __name__=="__main__":
    check_files(options.input)
