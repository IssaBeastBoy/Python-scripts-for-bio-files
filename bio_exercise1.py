
#print("Question 1\n\n")

with open("organism_1.fasta", "r") as fasta:
    fasta_content = fasta.read()
    fasta_content = fasta_content.split(">")
    fasta_content.remove("")
    for seq in fasta_content:
        seq = seq.split("\n")
        start = 1
        amino_acids = ""
        while(start < len(seq)):
            amino_acids = amino_acids + seq[start]
            start = start + 1
    header = seq[0]
    sequences = amino_acids

with open("report_1.txt","w+") as report:
    descri = header.split("|")
    report.write("Sequence descriptions:\n\t"+descri[4]+"\n")
    report.write("Sequence length:\n\tNucleotides - " + str(len(sequences)) + "\n")
    sequences = sequences.upper()
    nucleotides = list(sequences)
    ade = 0
    gua = 0
    thy = 0
    cyto = 0
    for nucleotide in nucleotides:
        if nucleotide == "A":
            ade = ade + 1
        elif nucleotide == "G":
            gua = gua + 1
        elif nucleotide == "T":
            thy = thy + 1
        elif nucleotide == "C":
            cyto = cyto + 1
        else:
            raise Exception("Not nucleotide")
    report.write("Nucleotides count:\n\tAdenine - " + str(ade) + "\n\t"+ "Guanine- " + str(gua) + "\n\t"+"Thymine- " + str(thy) + "\n\t"+ "Cytosine- " + str(cyto) + "\n\t")
    
#Question 2
def seq_count (header, sequences):
    with open("report_2.txt","a+") as report:
        descri = header.split("|")
        report.write("Sequence descriptions:\n\t"+descri[4]+"\n")
        report.write("Sequence length:\n\tNucleotides - " + str(len(sequences)) + "\n")
        sequences = sequences.upper()
        nucleotides = list(sequences)
        ade = 0
        gua = 0
        thy = 0
        cyto = 0
        for nucleotide in nucleotides:
            if nucleotide == "A":
                ade = ade + 1
            elif nucleotide == "G":
                gua = gua + 1
            elif nucleotide == "T":
                thy = thy + 1
            elif nucleotide == "C":
                cyto = cyto + 1
            else:
                raise Exception("Not nucleotide")
        report.write("Nucleotides count:\n\tAdenine - " + str(ade) + "\n\t"+ "Guanine- " + str(gua) + "\n\t"+"Thymine- " + str(thy) + "\n\t"+ "Cytosine- " + str(cyto) + "\n\n")    

with open("append_FASTA.fasta", "w+") as append:
    with open("organism_1.fasta", "r") as one:
        lines = one.readlines()
        for line in lines:
            append.write(line)
    with open("organism_2.fasta", "r") as two:
        lines = two.readlines()
        for line in lines:
            append.write(line)
    
with open("append_FASTA.fasta", "r") as  fasta:
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
for header, sequences in sequences.items():
    seq_count(header, sequences)


#Question 3
with open("hsp70_multi.fa", "r") as fasta:
    fasta_content = fasta.read()
    fasta_content = fasta_content.split(">")
    fasta_content.remove("")
    output = ""
    for seq in fasta_content:
        seq = seq.split("\n")
        output = output + seq[0]
        start = 1
        amino_acids = ""
        while(start < len(seq)):
            amino_acids = amino_acids + seq[start]
            start = start + 1
        output = output + ":\n\t Amino acid count: " +str(len(amino_acids)) +"\n"
    print("\nQuestion 3\n")
    print(output)
