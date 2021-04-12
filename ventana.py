from tkinter import  *
from tkinter import StringVar
from tkinter import filedialog
import os, sys, gzip, argparse, re, glob, numpy

root = Tk()
root.geometry("400x400")
root.title("MacroComplex Builder")


# def getFolderPath():
#     folder_selected = filedialog.askdirectory()
#     folderPath.set(folder_selected)
#
# def doStuff():
#     folder = folderPath.get()
#     print("Reading input files from", folder)

#
def run_program():
    verb = verbose_var.get()
    out_dir = output_dir
    max_file = max_files.get()
    # 1. Save all pdb files passed by the user in a list: files_list.
    files_list = self.file_name_paths
    print(files_list)
    folder = folderPath.get()
    print("Reading input files from", folder)

    if verb:
        sys.stderr.write("Program finished correctly. \n")

def getFolderPath():
    folder_selected = filedialog.askdirectory()
    os.chdir(folder_selected)
    dir = os.getcwd()
    dir2 = os.chdir("../")
    print(dir2)

folderPath = StringVar()
a = Label(root ,text="Input folder")
a.grid(row=0,column = 0)
E = Entry(root,textvariable=folderPath)
E.grid(row=0,column=1)
btnFind = Button(root, text="Browse Folder",command=getFolderPath)
btnFind.grid(row=0,column=2)

c = Button(root ,text="Run", command=run_program)
c.grid(row=4,column=0)

root.mainloop()



# from tkinter import  *
# from tkinter import StringVar
# from tkinter import filedialog
#
# root = Tk()
# root.geometry("400x400")
# root.title("MacroComplex Builder")
#
#

#
# def doStuff():
#     folder = folderPath.get()
#     os.path.isdir(folder)
#     # print("reading infp")
#     # #print(folder)
#     # for file in folder:
#     #     print(file)
#
# folderPath = StringVar()
# a = Label(root ,text="Input folder", font=("Arial Bold", 20))
# a.grid(row=0,column = 0)
# E = Entry(root,textvariable=folderPath)
# E.grid(row=0,column=1)
# btnFind = Button(root, text="Browse Folder",command=getFolderPath, font=("Arial Bold", 20))
# btnFind.grid(row=0,column=2)
#
# c = Button(root ,text="Run", command=doStuff)
# c.grid(row=4,column=0)
#
# root.mainloop()
