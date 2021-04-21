'''
with open("sequences.fa", "r") as fasta:
    fasta_content = fasta.read()
    fasta_content = fasta_content.split(">")
    fasta_content.remove("")
    sequences = {}
    for seq in fasta_content:
        seq = seq.split("\n")
        start = 1
        amino_acids = ""
        while(start < len(seq)):
            amino_acids = amino_acids + seq[start]
            start = start + 1
        sequences[seq[0]] = amino_acids
    count = 0
    losest = ["", 0,"", 0,"", 0 ]
    first_count= 0
    for key, value in sequences.items():
        if first_count == 0:
            losest[0] = key
            losest[1] =  len(value)
        if first_count == 1:
            losest[2] = key
            losest[3] =  len(value)
        if first_count == 2:
            losest[4] = key
            losest[5] = len(value)
        first_count = first_count + 1
        if first_count == 3:
            parse = 1
            next_seq = 1
            curr_key = losest[parse-1]
            curr_value = losest[parse]
            while True:
                if parse == 5:
                    parse = 1
                    curr_key = losest[parse-1]
                    curr_value = losest[parse]
                    next_seq = next_seq + 1
                if next_seq == 3:
                    break
                if curr_value > losest[parse + 2]:
                    losest[parse] = losest[parse + 2]
                    losest[parse-1] = losest[parse+1]
                
                    losest[parse + 2] = curr_value
                    losest[parse+1] = curr_key
                parse = parse + 2
        else:
            one = True
            two = True
            if  len(value) < losest[5]:                
                if len(value) < losest[3]:
                    temp_key = losest[2]
                    temp_value = losest[3]
                    losest[4] = temp_key
                    losest[5] = temp_value
                    if len(value) < losest[1]:
                        temp_key = losest[0]
                        temp_value = losest[1]
                        losest[2] = temp_key
                        losest[3] = temp_value
                        losest[0] = key
                        losest[1] =  len(value)
                        one = False
                    if one:
                        losest[2] = key
                        losest[3] =  len(value)
                        two = False
                if two and one:
                    losest[4] = key
                    losest[5] =  len(value)
print("Question 1\n")
print(losest)

with open("sequences.fa", "r") as fasta:
    fasta_content = fasta.read()
    fasta_content = fasta_content.split(">")
    fasta_content.remove("")
    sequences = {}
    for seq in fasta_content:
        seq = seq.split("\n")
        start = 1
        amino_acids = ""
        while(start < len(seq)):
            amino_acids = amino_acids + seq[start]
            start = start + 1
        sequences[seq[0]] = amino_acids
    output = ""
    for key, value in sequences.items():
        output = output + key + "\n"
        length = len(value)
        start = 0
        end = 60        
        while (end > 0):
            output = output + value[start:end] + "\n"
            length = length - 60
            if 60%length >= 60:
                start = end
                end = end + 60
            else:
                start = end
                end = len(value)
                output = output + value[start:end]
                end = -1
        output = output + "\n\n"
        
with open("sequences_60.fa", "w") as fasta:
    fasta.write(output)

with open("alignment_Hsp70.txt") as converter: 
    sequences = converter.readlines()
    output_list = []
    first_parse = 0
    header = 0
    for line in sequences:
        probe = list(line)
        if probe[0] != ' ' and probe[0] != '\n' and first_parse == 0:
            sequence = line.split(" ")
            store = ""
            store = store + ">" + sequence[0]
            if(sequence[1] !=''):
                start = 2
                store = store + ' ' + sequence[1]
            else:
                start = 1                       
            while(start < len(sequence)):
                if len(sequence[start]) > 6:
                    store = store + "\n" + sequence[start]
                    break
                start = start+1
            output_list.append(store)
        
        if probe[0] != ' ' and probe[0] != '\n' and first_parse > 1:
            sequence = line.split(" ")   
            store = ""         
            start = 1            
            while(start < len(sequence)):
                if len(sequence[start]) > 15:
                    store = "\n" + sequence[start]
                    break
                start = start + 1
            temp = output_list[header]
            store = temp + store + "\n"
            output_list[header] = store
            header = header + 1
            
        if probe[0] == '\n':
             first_parse = first_parse + 1
             header = 0
    with open("alignment_Hsp70.fa", "w") as aligned:
        for sequence in output_list:
            aligned.write(sequence)
 '''

with open("alignment_Hsp70.txt") as converter: 
    sequences = converter.readlines()
    output_list = []
    first_parse = 0
    header = 0
    for line in sequences:
        probe = list(line)        
        if probe[0] != ' ' and probe[0] != '\n' and first_parse == 0:
            line = ""
            for chars in probe:
                if chars != '-':
                    line = line + chars
            sequence = line.split(" ")
            store = ""
            store = store + ">" + sequence[0]
            if(sequence[1] !=''):
                start = 2
                store = store + ' ' + sequence[1]
            else:
                start = 1                       
            while(start < len(sequence)):
                if len(sequence[start]) > 6:
                    store = store + "\n" + sequence[start]
                    break
                start = start+1
            output_list.append(store)
        
        if probe[0] != ' ' and probe[0] != '\n' and first_parse > 1:
            line = ""
            for chars in probe:
                if chars != '-':
                    line = line + chars
            sequence = line.split(" ")   
            store = ""         
            start = 1            
            while(start < len(sequence)):
                if len(sequence[start]) > 15:
                    store = "\n" + sequence[start]
                    break
                start = start + 1
            temp = output_list[header]
            store = temp + store + "\n"
            output_list[header] = store
            header = header + 1
            
        if probe[0] == '\n':
             first_parse = first_parse + 1
             header = 0
    with open("Hsp70.fa", "w") as aligned:
        for sequence in output_list:
            aligned.write(sequence)      
                
            

        

        