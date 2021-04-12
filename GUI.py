from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import Checkbutton
from tkinter import StringVar
from tkinter import Entry
import os, sys, gzip, argparse, re, glob, numpy
from functions import *
from Bio.PDB import *


class Window(Frame):

    def client_exit(self):
        if messagebox.askyesno("EXIT", "Are you sure you want to quit the program?"):
            Frame.quit(self)

    def choosedir(self):
        self.folder_path = filedialog.askdirectory(title="Select folder with desired PDB files")
        if self.folder_path:
            for file in os.listdir(self.folder_path):
                name = re.search('.*\.pdb', file)
                if name:
                    file_name = self.folder_path + "/" + name.group()
                    self.file_name_paths.append(file_name)
                    self.list_files.append(name.group())

        self.list_filesbox.delete(0, END)
        index = 1
        for record in self.list_files:
            self.list_filesbox.insert(index, record)
            index += 1

    def select_outdir(self):
        self.output_dir = filedialog.askdirectory(title="Select desired output directory")

    def create_left_frame(self):
        left_frame = LabelFrame(self, text="Files List", width=150, padx=5, pady=5)

        frame = Frame(left_frame)

        scrollbar = Scrollbar(frame, orient=VERTICAL)

        self.list_filesbox = Listbox(frame, selectmode=SINGLE, height=20, yscrollcommand=scrollbar.set)

        scrollbar.config(command=self.list_filesbox.yview)
        scrollbar.pack(side=RIGHT, fill=Y)

        self.list_filesbox.pack(side=LEFT, expand=True, fill=BOTH)

        frame.pack(fill=BOTH)

        run_frame = Frame(left_frame, width=150)

        b= Button(run_frame, text="Run program!", command=self.runProgram)
        b.pack(side=RIGHT)
        b2 = Button(run_frame, text="Clean", command=self.clean)
        b2.pack()

        run_frame.pack(fill=BOTH)

        left_frame.grid(row=0, column=0)


    def create_output_frame(self):
        self.output_frame = LabelFrame(self, text="Output Files", width=400, padx=5, pady=5)

        output_options_frame = Frame(self.output_frame, width=390, borderwidth=2)

        self.verbose_var = BooleanVar(value=False)
        self.st_var = BooleanVar(value=False)

        b2 = Checkbutton(output_options_frame, text="Verbose", variable=self.verbose_var, onvalue=True, offvalue=False)
        b2.pack(side=LEFT)

        output_options_frame.pack(fill=BOTH)

        second_options_frame = Frame(self.output_frame, width=390, borderwidth=2)

        self.max_files = IntVar(value=10)

        L1 = Label(second_options_frame, text="Enter maximum number files")
        L1.pack(side=LEFT)
        max_files_entry = Entry(second_options_frame, width=5, textvariable=self.max_files)
        max_files_entry.pack(side=RIGHT)

        second_options_frame.pack(fill=BOTH)

        third_options_frame = Frame(self.output_frame, width=390, borderwidth=2)

        L2 = Label(third_options_frame, text="Enter Desired output directory")
        L2.pack(side=LEFT)
        output_dir_entry = Button(third_options_frame, text="Browse", command=self.select_outdir)
        output_dir_entry.pack(side=RIGHT)

        third_options_frame.pack(fill=BOTH)

        fourth_options_frame = Frame(self.output_frame, width=390, borderwidth=2)

        self.st_var = StringVar(value="A1B1")
        L3 = Label(fourth_options_frame, text="Enter Stoichiometry")
        L3.pack(side=LEFT)
        st_entry = Entry(fourth_options_frame, width=20, textvariable=self.st_var)
        st_entry.pack(side=RIGHT)

        fourth_options_frame.pack(fill=BOTH)

        b4 = Button(self.output_frame, text="Show options", command=self.showOptions)
        b4.pack()

        #self.output_complex = Text(self.output_frame)
        #self.output_complex.pack()

        self.output_frame.grid(row=0, column=1)

    def showOptions(self):

        self.message = StringVar()
        self.message.set("Verbose: %s, Stoichiometry: %s, output directory: %s, maximum files: %s" % (
            self.verbose_var.get(), self.st_var.get(), self.output_dir, self.max_files.get()))
        messagebox.showinfo("checkbox values", self.message.get())


    def clean(self):
        self.output_dir.set("")
        self.max_files.set("")
        self.st_var.set(0)
        self.verbose_var.set(0)


    def save_final_PDB(self):
        self.output_name = filedialog.asksaveasfilename(title="Save Final complex as", defaultextension=".pdb")

    def choosefiles(self):
        st_file = filedialog.askopenfilename(title="Select file with Stoichiometry", filetypes=(("Text files", "*.txt"),("all files", "*.*")))
        self.st_filepath = st_file

    def showinfo(self):
        messagebox.showinfo(title="Info", message="This program has been done for SBI and PYT subject projects. Authors: Ramon Massoni and Winona Oliveros")

    def showhelp(self):
        messagebox.showinfo(title="Help", message="To see the documentation of this app and the tutorial go to: \n https://github.com/massonix/Python-project")

    def runProgram(self):

        verb = self.verbose_var.get()
        out_dir = self.output_dir
        max_file = self.max_files.get()

        # 1. Save all pdb files passed by the user in a list: files_list.

        files_list = self.file_name_paths

        # 2. Create a  dictionary with the files path as key and the structure object as values: structures.

        structures = {}
        ids = list(range(1000))
        i = 0
        for f in files_list:
            structures[f] = GetStructures(f, ids[i:i + 2])
            i += 2

        if verb:
            sys.stderr.write("%d PDB files found: \n" % len(structures))
            for f in structures.keys():
                sys.stderr.write(f + "\n")

        # 3. Process the stoichiometry input to calculate the recurssion depth (k) and create a dictionary: stoich_dict.

        stoich_in = self.st_var.get()
        stoich_dict = {}
        i = 0
        j = 0
        while i < len(stoich_in):
            if stoich_in[i].isalpha():
                stoich_dict[stoich_in[i]] = ""
                j += 1
                while stoich_in[j].isdigit():
                    stoich_dict[stoich_in[i]] += str(stoich_in[j])
                    if j < (len(stoich_in) - 1):
                        j += 1
                    else:
                        break
            i += 1
            j = i

        stoich_dict = {x: int(y) for x, y in stoich_dict.items()}
        k = sum(list(stoich_dict.values()))

        if verb:
            sys.stderr.write("Stoichiometry: %s\n" % stoich_in)

        # 4. Map all the similar chains (>95% sequence similarity) in a dictionary: similar_chains.
        # The keys are the chain ids, and the values the id of the first chain they are similar to.

        structures2 = structures.copy()
        scores = {}

        structures_iter = list(structures2.items())
        similar_chains = {}
        i = 0

        for index1, items1 in enumerate(structures_iter):
            file1 = items1[0]
            structure1 = items1[1]
            chains1 = list(structure1[0].get_chains())
            for file2, structure2 in structures_iter[index1:len(structures_iter)]:
                chains2 = list(structure2[0].get_chains())

                for chain1 in chains1:
                    i += 1
                    for chain2 in chains2:
                        if chain2.id in similar_chains:
                            continue
                        i += 1
                        Alignment = Alignsequence(chain1, chain2)
                        score = Alignment[0][2] / len(Alignment[0][0])
                        if score > 0.95:
                            similar_chains.setdefault(chain2.id, chain1.id)

        # 5. Remove those structures that do not share any similar chain with another structure and thus cannot be superimposed.

        similar_chains_keys = sorted(list(similar_chains.keys()))
        similar_chains_values = [similar_chains[x] for x in similar_chains_keys]
        similar_chains_keys_iter = [(similar_chains_keys[x], similar_chains_keys[x + 1]) for x in
                                    range(0, len(similar_chains_keys), 2)]
        similar_chains_values_iter = [(similar_chains_values[x], similar_chains_values[x + 1]) for x in
                                      range(0, len(similar_chains_values), 2)]

        for i, chs in enumerate(similar_chains_values_iter):
            if similar_chains_values.count(chs[0]) == 1 and similar_chains_values.count(chs[1]) == 1:
                structures.pop(files_list[i])
                del similar_chains[similar_chains_keys_iter[i][0]]
                del similar_chains[similar_chains_keys_iter[i][1]]
            elif chs[0] == chs[1] and similar_chains_values.count(chs[0]) == 2:
                structures.pop(files_list[i])
                del similar_chains[similar_chains_keys_iter[i][0]]
                del similar_chains[similar_chains_keys_iter[i][1]]

        if not similar_chains:
            raise ValueError("Unable to superimpose: no common chains among structures.")

        del structures2
        del scores
        del structures_iter
        del similar_chains_keys
        del similar_chains_values
        del similar_chains_keys_iter

        # 6. Check if there are enough different chains to achieve the desired stoichiometry.

        stoich_set = set(list(similar_chains.values()))

        if len(stoich_set) < len(stoich_dict):
            raise ValueError("Impossibe stoichiometry: The provided stoichiometry contains %d different chains, but "
                             "the input PDBs only have %d unique chains." % (len(stoich_dict), len(stoich_set)))

        # 7. Call the recursive function and build the complexes with the desired stoichiometry.

        contact_num = False

        if contact_num:
            contacts = contact_num
        else:
            contacts = 5

        strs = [structures[x] for x in files_list]
        stoich_list = list(stoich_set)
        strs_similar_dic = {strs[x]: similar_chains_values_iter[x] for x in range(len(strs))}

        for ind, uniq_chain1 in enumerate(stoich_list):
            for uniq_chain2 in stoich_list[ind:]:
                dimer = [x for x in strs if
                         strs_similar_dic[x] == (uniq_chain1, uniq_chain2) or strs_similar_dic[x] == (
                             uniq_chain2, uniq_chain1)]
                if dimer:
                    if uniq_chain1 == uniq_chain2:
                        if len(stoich_dict) == 1:
                            strs2 = [x for x in strs if strs_similar_dic[x] == (uniq_chain1, uniq_chain2)]
                            if verb:
                                sys.stderr.write("Starting recursion...\n")
                                print(strs2[0].get_chains())
                            impose_clash(dimer[0], strs2, k, 2, contacts, 0, similar_chains, stoich_dict, out_dir,
                                         max_file, verb)
                        else:
                            continue
                    else:
                        if len(stoich_dict) == 1:
                            continue
                        else:
                            if verb:
                                sys.stderr.write("Starting recursion...\n")
                            impose_clash(dimer[0], strs, k, 2, contacts, 0, similar_chains, stoich_dict, out_dir,
                                         max_file, verb)


        if verb:
            sys.stderr.write("Program finished correctly. \n")

    def create_menu(self):
        self.menu = Menu(self)

        #FOLDER MENU
        folder_menu = Menu(self.menu)
        folder_menu.add_command(label="Select Folder", command=self.choosedir)
        folder_menu.add_command(label="Select Stoichiometry file", command=self.choosefiles)
        folder_menu.add_separator()
        folder_menu.add_command(label="Save Final PDB as", command=self.save_final_PDB)
        folder_menu.add_separator()
        folder_menu.add_command(label="Exit", command=self.client_exit)

        #HELP MENU
        help_menu = Menu(self.menu)
        help_menu.add_command(label="Information", command=self.showinfo)
        help_menu.add_command(label="Help", command=self.showhelp)

        self.menu.add_cascade(label="Folder/File Opts", menu=folder_menu)
        self.menu.add_cascade(label="Help", menu=help_menu)

        self.master.config(menu=self.menu)

    def createWidgets(self):
        self.create_menu()
        self.create_left_frame()
        self.create_output_frame()

        self.grid(row=0)

    def __init__(self, master=None, **kwargs):
        Frame.__init__(self, master, **kwargs)

        self.master = master
        self.master.title("Assemblator plus")
        # Impide que los bordes puedan desplazarse para
        # ampliar o reducir el tamaño de la ventana 'self.raiz':
        self.master.resizable(width=False, height=False)

        self.config(width=600)
        self.config(width=600)

        self.menu = None
        self.list_files = []
        self.folder_path = None
        self.list_filesbox = None
        self.output_complex = None
        self.st_filepath = None
        self.output_name = None
        self.left_frame = None
        self.file_name_paths = []
        self.max_files = None
        self.st_var = None
        self.output_dir = None
        self.message = None
        self.verbose_var = None

        self.createWidgets()


root = Tk()


app = Window(master=root, padx=10, pady=10)

app.mainloop()
