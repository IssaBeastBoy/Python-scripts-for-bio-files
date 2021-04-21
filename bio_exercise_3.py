
with open("1yuw.pdb", "r") as structure:
    protein_Info = {}
    contents = structure.readlines()
    skip = 0
    amino_sequence = ""
    atom_information = ""
    for probe in contents:
        line = probe.split(' ')
        if line[0] == "SEQRES":
            one_AMINO = ["R", "K", "D", "E", "Q", "N", "H", "S", "T", "Y", "C", "W", "M", "A", "I", "L", "F", "V", "P", "G"]
            three_AMINO = ["ARG", "LYS", "ASP", "GLU", "GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS", "TRP", "MET", "ALA", "ILA", "LEU", "PHE", "VAL", "PRO", "GLY"]
            for amino in line:
                amino = amino.split("\n")
                if three_AMINO.__contains__(amino[0]):
                    index = three_AMINO.index(amino[0])
                    amino_sequence = amino_sequence + one_AMINO[index]
        elif line[0] == "ATOM":
            start = 1
            while start < len(line):
                if line[start] != '':
                    atom_information = atom_information + " " + line[start]
                start = start + 1
        elif line[0] == "TER":
            start = 1
            while start < len(line):
                if line[start] != '':
                    atom_information = atom_information + " " + line[start]
                start = start + 1
        else:
            start = 0
            if skip < 4:
                for parts in line:
                    if parts == "TITLE":
                        found = line[start:len(probe)]
                        store = ""
                        start = 1
                        
                        while start < len(found):
                            if found[start] != '':
                                break
                            start = start + 1
                        
                        while start < len(found):
                            store = store + " " + found[start]
                            start = start + 1
                        protein_Info["Title"] = store
                        skip = skip + 1
                        break
                    elif parts == "MOL_ID:":
                        found = line[start:len(probe)]
                        store = ""
                        for output in found:
                            store = store + output
                        protein_Info["ID"] = store
                        break
                    elif parts == "MOLECULE:":
                        found = line[start:len(probe)]
                        store = ""
                        for output in found:
                            if output == '':
                                store = store + "\n"
                                break
                            else:
                                store = store + output + " "
                        protein_Info["Associated Mol"] = store
                        skip = skip + 1
                        break
                    elif parts == "SYNONYM:":
                        found = line[start:len(probe)]
                        store = ""
                        start = 1
                        while start < len(found):
                            store = store + " " + found[start]
                            start = start + 1
                        protein_Info["Synonyms"] = store
                        skip = skip + 1
                        break
                    else:
                        start = start + 1
    protein_Info["Protein sequence"] = amino_sequence
    protein_Info["Atom"] = atom_information


with open("4hhb.pdb", "r") as structure:
    protein_Info = {}
    contents = structure.readlines()
    skip = 0
    amino_sequence = ""
    atom_information = ""
    title = 0
    mol = 0
    mol_list = []
    chains = []
    for probe in contents:
        line = probe.split(' ')
        start = 0
        if skip == 0:
            for parts in line:
                if parts == "TITLE" and title == 0:
                    found = line[start:len(probe)]
                    store = ""
                    start = 1
                    title = title + 1
                    
                    while start < len(found):
                        if found[start] != '':
                            break
                        start = start + 1
                    
                    while start < len(found):
                        store = store + " " + found[start]
                        start = start + 1
                    protein_Info["Title"] = store
                    break
                elif parts == "MOL_ID:":
                    found = line[start:len(probe)]
                    found = found[0:2]
                    store = ""
                    for output in found:
                        store = store + output
                    if protein_Info.get(store) == None:
                        mol = mol + 1
                        mol_list.append(store)
                        protein_Info[store] = {}
                    break
                elif parts == "MOLECULE:":
                    found = line[start:len(probe)]
                    store = ""
                    start = 1
                    key = "MOL_ID:" + str(mol) +";"
                    while start < len(found):
                        if found[start] == '':
                            store = store + "\n"
                            break
                        else:
                            store = store + found[start] + " "
                        start = start + 1
                    protein_Info[key]["Associated Mol"] = store
                    break
                elif parts == "SOURCE":
                    skip = skip + 1
                    break
                elif parts == "CHAIN:":
                    found = line[start:len(probe)]
                    store = ""
                    start = 1
                    key = "MOL_ID:" + str(mol) +";"
                    while start < len(found):
                        if found[start] =='':
                            break
                        chain = list(found[start])
                        protein_Info[key][chain[0]] = ""
                        chains.append(chain[0])
                        start = start + 1
                    break
                else:
                    start = start + 1
        else:
            for seq in line:
                if seq == "SEQRES":
                    one_AMINO = ["R", "K", "D", "E", "Q", "N", "H", "S", "T", "Y", "C", "W", "M", "A", "I", "L", "F", "V", "P", "G"]
                    three_AMINO = ["ARG", "LYS", "ASP", "GLU", "GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS", "TRP", "MET", "ALA", "ILA", "LEU", "PHE", "VAL", "PRO", "GLY"]            
                    for mol in mol_list:
                        information = protein_Info.get(mol)
                        start = 1
                        while start < len(line):
                            if chains.__contains__(line[start]):
                                break
                            start = start + 1
                        if information.get(line[start]) != None:
                            amino_parse = start
                            sequence = ""
                            while amino_parse < len(line):
                                if three_AMINO.__contains__(line[amino_parse]):
                                    index = one_AMINO[three_AMINO.index(line[amino_parse])]
                                    sequence = sequence + index
                                amino_parse = amino_parse + 1
                            temp = protein_Info[mol][line[start]]
                            temp = temp + sequence
                            protein_Info[mol][line[start]] = temp                            
                            break
                elif seq == "ATOM":
                    start = 1
                    while start < len(line):
                        if line[start] != '':
                            atom_information = atom_information + " " + line[start]
                        start = start + 1
                elif seq == "TER":
                    start = 1
                    while start < len(line):
                        if line[start] != '':
                            atom_information = atom_information + " " + line[start]
                        start = start + 1
                else:
                    break
    protein_Info["Atom"] = atom_information
    


