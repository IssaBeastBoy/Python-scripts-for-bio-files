import os
import re

def read_AMINO(content, header, type, atom):
    sequence = ""  
    backbone = ""
    alphaCO2 = ""
    output = []
    if type == 0: 
        for lines in content:
            line = lines.split(" ")
            start = 1  
            while start < len(line):
                if line[0].__contains__(header):
                    if line[start] != "" and len(line[start]) > 2:
                        try:
                            amino = line[start].split("\n")
                            number = int(amino[0])
                        except:
                            amino = line[start].split("\n")
                            if atom == -1:
                                one_AMINO = ["R", "K", "D", "E", "Q", "N", "H", "S", "T", "Y", "C", "W", "M", "A", "I", "L", "F", "V", "P", "G"]
                                three_AMINO = ["ARG", "LYS", "ASP", "GLU", "GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS", "TRP", "MET", "ALA", "ILA", "LEU", "PHE", "VAL", "PRO", "GLY"]
                                if three_AMINO.__contains__(amino[0]):
                                    sequence = sequence + one_AMINO[three_AMINO.index(amino[0])]
                            else:
                                sequence = sequence + amino[0]   
                start = start + 1          
        return sequence
    if type == 1:
        for lines in content:
            line = lines.split(" ")
            start = 1 
            checker1 = 0
            checker2 = 0
            temp = ""
            while start < len(line):
                if line[0].__contains__(header):
                    if atom == -3:
                        one_AMINO = ["R", "K", "D", "E", "Q", "N", "H", "S", "T", "Y", "C", "W", "M", "A", "I", "L", "F", "V", "P", "G"]
                        three_AMINO = ["ARG", "LYS", "ASP", "GLU", "GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS", "TRP", "MET", "ALA", "ILA", "LEU", "PHE", "VAL", "PRO", "GLY"]
                        if three_AMINO.__contains__(line[start]):
                            if len(output)> 0:
                                if line[start] != three_AMINO[one_AMINO.index(output[checker1 - 1])]:
                                    output.append(one_AMINO[three_AMINO.index(line[start])])
                                    checker1 = checker1 + 1
                            else:
                                output.append(one_AMINO[three_AMINO.index(line[start])])
                                checker1 = checker1 + 1
                    else:
                        if (line[start] == "CA" or line[start] == "C" or line[start] == "N" or line[start] == "O") and start <len(line) - 5:
                            checker1 = 1                
                        if (line[start] == "CA"):
                            checker2 = 1
                        if line[start] == "":
                            temp = temp + " " 
                        else:
                            temp = temp + " " + line[start]
                start = start + 1             
            if line[0] == header:
                if checker1 == 1:
                    backbone = backbone + temp + "\n"
                if checker2 == 1:
                    alphaCO2 = alphaCO2 + temp + "\n"
                sequence = sequence + temp  + "\n"
    if atom == 0:
        return backbone
    elif atom == 1:
        return alphaCO2
    elif atom == -3:
        return output
    else:
        return sequence

def info_position(matches, type, sequence):
    positions = []
    for probe in matches:
        positions.append([probe.start(), probe.end()])
    pos_output = ""
    if positions == []:
        return type + " sequence not found"
    for pos in positions:
        motif_pos = sequence[pos[0]:pos[1]]
        pos_output = pos_output + "\n\tstart position - "+ str(pos[0]) + " end position - " + str(pos[1])
        motif_pos = motif_pos.lower()
        store = []
        if pos[0] > 0:
            index = 0
            store.append(sequence[index: pos[0]])
            index = len(sequence)
            store.append(sequence[pos[1]: index])
            sequence = store[0]+motif_pos+store[1] 
        else:
            sequence = motif_pos + sequence[pos[1]:len(sequence)]
    output = type + " appears " + str(len(positions)) + pos_output + "\n" + sequence
    return output

read_files = {}
print("Welcome to PDE process program please select the given options: ")
file_name = ''
while True:
    user_input = input("(R) - Read \n(S) - Search \n(W) - Write \n(I) - Information \n(H) - Help \n(Q) - Quit \n")
    breaker = 0
    if user_input == "Q":
        break
    elif  user_input == "I":
        if read_files:
            option = input("\nDisplays FASTA sequences (Y/N): ")
            if option == "Y":
                print("\n\tChoose the given files blow:")
                print()
                for key in read_files:
                    print("\t*"+key)
                name = input("\n\tEnter file name: ")
                print()
                if read_files.get(name) != None:
                    option = input("\n\t\tDisplay SEQRES sequence (Y/N): ")
                    if option == "Y":
                            title = name + "_SEQRES.fa"
                            content = read_files.get(name)
                            content = content.split("\n")
                            sequence = read_AMINO(content, "SEQRES", 0,-1)
                            sequence = list(sequence)
                            start = 0
                            output = "> " + name +"\n"
                            for amino in sequence:
                                if start == 61:
                                    start = 0
                                    output = output + "\n"
                                output = output + amino
                                start = start + 1
                            print(output)
                            print("\n")
                    option = input("\n\t\tCoordinates sequence FASTA File (Y/N): ")
                    if option == "Y":
                        title = name + "_Coordinates.fa"
                        content = read_files.get(name)
                        content = content.split("\n")
                        sequence = read_AMINO(content, "ATOM", 1,-3)
                        start = 0
                        output = "> " + name +"\n"
                        for amino in sequence:
                            if start == 61:
                                start = 0
                                output = output + "\n"
                            output = output + amino
                            start = start + 1
                        print(output)
                        print("\n")
                    option = input("\n\t\tAligment sequence FASTA File (Y/N): ")
                    if option == "Y":
                        title = name + "_Alignment.fa"
                        content = read_files.get(name)
                        content = content.split("\n")
                        ATOM = read_AMINO(content, "ATOM", 1,-3)
                        SEQRES = read_AMINO(content, "SEQRES", 0,-1)
                        SEQRES = list(SEQRES)
                        start = 0
                        index = 0
                        output = "> " + name +"\n"
                        store = ""
                        while(index < len(ATOM)):
                            if start < len(SEQRES):
                                if ATOM[index] == SEQRES[start]:
                                    store = store + ATOM[index]
                                    index = index + 1
                                    start = start + 1
                                else:
                                    store = store + "X"
                                    start = start + 1
                            else:
                                store = store + ATOM[index]
                                index = index + 1
                        sequence = list(store)
                        start = 0
                        for amino in sequence:
                            if start == 61:
                                start = 0
                                output = output + "\n"
                            output = output + amino
                            start = start + 1
                        print(output)
                        print("\n")
                    option = input("\n\t\tDisplay non-water ligands sequence (Y/N): ")
                    if option == "Y":
                        content = read_files.get(name)
                        content = content.split("\n")
                        output = ""
                        for lines in content:
                            line = lines.split(" ")
                            if line[0].__contains__("HETATM"):
                                if not line.__contains__("HOH"):
                                    store = ""
                                    for chars in line:
                                        if chars == "":
                                            store = store + " "
                                        else:
                                            store = store + " " + chars 
                                output = output + store +"\n"
                        print(output)
                        print("\n")
                else:
                    print("No such file in storeage")
        else:
            print("\n Error - No files read in \n")                
    elif user_input == "H":
        output = "\n\nR - Read in file in folder\
                \nS - Search charactertics files which are already read in\
                    \n\t* Motif search: Search for motif sites in selected file\
                    \n\t* Glycosylation search: Search for glycosylation sites in selected file\
                \nW - Write selected information into a new file\
                    \n\t Write Atoms - Write information regarding atoms to a new file\
                        \n\t\t* All atom: All atom coordinates to new file\
                        \n\t\t* Backbone atoms: All backbone atoms coordinates to new file\
                        \n\t\t* Alpha carbon atoms: All alpha carbon atoms coordinates to new file\
                    \n\t Write Fasta - Write sequence in fasta format\
                        \n\t\t* SEQRES Sequence: All SEQRES sequence to new file\
                        \n\t\t* Coordinate Sequence: All coordinate sequence to new file\
                        \n\t\t* Alignment Sequence: All alignment sequence to new file\
                \nI - Information being print out onto the console\
                    \n\t* SEQRES Sequence: Print SEQRES sequence to console\
                    \n\t* Coordinate Sequence: Print Ccoordinate sequence to console\
                    \n\t* Alignment Sequence: Print alignment sequence to console\
                    \n\t* Non water ligands: Print hydrophilic ligands\
                \nH - Help will give information regards all functionality offered\
                \nQ - Quite program\n\n"
        print(output)
    elif user_input == "R":
        file_name = input("Please enter file name: ")
        print()
        path = os.getcwd()
        os.chdir(path)
        folder = os.listdir()
        parse = 0
        while parse < len(folder):
            if file_name == folder[parse]:
                if read_files.get(file_name) == None:
                    with open(file_name) as content:
                        information = content.read()
                        read_files[file_name] = information
                    breaker = 1
                else:
                    breaker = 1
                    print("\nFile already exist \n")
            parse = parse + 1
            if parse == len(folder) and breaker == 0:
                print("\nFile name does not exist\n")
                break
    elif user_input == "S":
        if file_name == '':
            print("\nNo files to select\n")
        else:
            print("\nChoose the given files blow:")
            print()
            for key in read_files:
                print(key)
            name = input("\nEnter file name: ")
            print()
            if read_files.get(name) != None:
                respond = input("\tCheck for motif(Y/N): ")
                if respond == "Y":
                    motif = input("\n\tMotif search - Enter sequence that you are looking for (e.g PRO GLN VAL THR LEU ): ")
                    checker = 1
                    parse = 0
                    temp = ""
                    tempM = motif.split(" ")
                    three_AMINO = ["ARG", "LYS", "ASP", "GLU", "GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS", "TRP", "MET", "ALA", "ILA", "LEU", "PHE", "VAL", "PRO", "GLY"]
                    for protein in tempM:
                        if three_AMINO.__contains__(protein):
                            temp = temp + protein
                        else:
                            print("\nERROR - incorrect amino acid\n")
                            checker = 0
                            break
                    if checker == 1:
                        motif = temp
                        content = read_files.get(name)
                        content = content.split("\n")
                        sequence = read_AMINO(content, "SEQRES", 0, -2)                        
                        matches = re.finditer(motif, sequence)
                        print (info_position(matches,"Motif", sequence))
                print()
                respond = input("\tCheck for glycosylation sites(Y/N): ")
                if respond == "Y":
                    glycosylation = "(G|Y|L)+A(P|F|W)+(T|L|V|M|I)"
                    content = read_files.get(name)
                    content = content.split("\n")
                    sequence = read_AMINO(content, "SEQRES", 0, -1)   
                    print("\n")
                    matches = re.finditer(glycosylation, sequence)
                    print (info_position(matches,"Glycoslyation", sequence))
            else:
                print("Incorrect File selected")
    elif user_input == "W":
        if read_files:
            option = input("\tGenerate coordinate File (Y/N): ")
            if option == "Y":            
                print("\n\tChoose the given files blow:")
                print()
                for key in read_files:
                    print("\t*"+key)
                name = input("\n\tEnter file name: ")
                print()
                if read_files.get(name) != None:
                    option = input("\n\t\tWrite all atom (Y/N): ")
                    if option == "Y":
                        title = name + "_All_atom.txt"
                        with open(title,"+w") as atoms:
                            content = read_files.get(name)
                            content = content.split("\n")
                            atom = read_AMINO(content, "ATOM", 1,-1)
                            atoms.write(atom)
                    option = input("\n\t\tWrite backbone atoms (Y/N): ")
                    if option == "Y":
                        title = name + "_backbone_atom.txt"
                        with open(title,"+w") as atoms:
                            content = read_files.get(name)
                            content = content.split("\n")
                            atom = read_AMINO(content, "ATOM", 1, 0)
                            atoms.write(atom)
                    option = input("\n\t\tWrite alpha carbon atoms (Y/N): ")
                    if option == "Y":
                        title = name + "_alphaCarbon_atom.txt"
                        with open(title,"+w") as atoms:
                            content = read_files.get(name)
                            content = content.split("\n")
                            atom = read_AMINO(content, "ATOM", 1, 1)
                            atoms.write(atom)
                
            option = input("\n\tGenerate FASTA File (Y/N): ")
            if option == "Y":
                print("\n\tChoose the given files blow:")
                print()
                for key in read_files:
                    print("\t*"+key)
                name = input("\n\tEnter file name: ")
                print()
                if read_files.get(name) != None:
                    option = input("\n\t\tSEQRES sequence FASTA File (Y/N): ")
                    if option == "Y":
                            title = name + "_SEQRES.fa"
                            with open(title,"+w") as S_fasta:
                                content = read_files.get(name)
                                content = content.split("\n")
                                sequence = read_AMINO(content, "SEQRES", 0,-1)
                                sequence = list(sequence)
                                start = 0
                                output = "> " + name +"\n"
                                for amino in sequence:
                                    if start == 61:
                                        start = 0
                                        output = output + "\n"
                                    output = output + amino
                                    start = start + 1
                                S_fasta.write(output)
                    option = input("\n\t\tCoordinates sequence FASTA File (Y/N): ")
                    if option == "Y":
                        title = name + "_Coordinates.fa"
                        with open(title,"+w") as C_fasta:
                            content = read_files.get(name)
                            content = content.split("\n")
                            sequence = read_AMINO(content, "ATOM", 1,-3)
                            start = 0
                            output = "> " + name +"\n"
                            for amino in sequence:
                                if start == 61:
                                    start = 0
                                    output = output + "\n"
                                output = output + amino
                                start = start + 1
                            C_fasta.write(output)
                    option = input("\n\t\tAligment sequence FASTA File (Y/N): ")
                    if option == "Y":
                        title = name + "_Alignment.fa"
                        with open(title,"+w") as A_fasta:
                            content = read_files.get(name)
                            content = content.split("\n")
                            ATOM = read_AMINO(content, "ATOM", 1,-3)
                            SEQRES = read_AMINO(content, "SEQRES", 0,-1)
                            SEQRES = list(SEQRES)
                            start = 0
                            index = 0
                            output = "> " + name +"\n"
                            store = ""
                            while(index < len(ATOM)):
                                if start < len(SEQRES):
                                    if ATOM[index] == SEQRES[start]:
                                        store = store + ATOM[index]
                                        index = index + 1
                                        start = start + 1
                                    else:
                                        store = store + "X"
                                        start = start + 1
                                else:
                                    store = store + ATOM[index]
                                    index = index + 1
                            sequence = list(store)
                            start = 0
                            for amino in sequence:
                                if start == 61:
                                    start = 0
                                    output = output + "\n"
                                output = output + amino
                                start = start + 1
                            A_fasta.write(output)
                else:
                    print("No such file in storeage")                
        else:
            print("\n Error - No files read in \n")
    