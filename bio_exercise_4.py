import re

with open("Pf3D7_08_v3_nt.fa","r") as gen_seq:
    gene = gen_seq.read()
    seq = gene.split("\n")
    start = 1
    chromo = ""
    while start < len(seq):
        chromo = chromo + seq[start]
        start = start + 1
    with open("PfHsp70-x_nt.fa","r") as hsp70_seq:
        gene = hsp70_seq.read()
        seq = gene.split("\n")
        start = 1
        sequence = ""
        mRNA_sequence = ""
        while start < len(seq):
            nucleotides = list(seq[start])
            for nucleotide in nucleotides:
                mRNA_sequence = mRNA_sequence + nucleotide
                if nucleotide =="A":
                    sequence = "T" + sequence
                elif nucleotide == "T":
                    sequence = "A" + sequence
                elif nucleotide == "C":
                    sequence = "G" + sequence
                elif nucleotide == "G":
                    sequence = "C" + sequence
                else:
                    raise Exception("Not nucleotide character")
            start = start + 1
    search = re.search(sequence,chromo)
    if(search):
        print("\nNucleotide sequence found")
        print("\nStart position: ", search.start(), "\n")
        print("End position: ", search.end(), "\n")
    else:
        print("Not found")
    expression = "(ATG[TCGA]+TAG)|(ATG[TCGA]+TAA)|(ATG[TCGA]+TGA)"
    object = re.compile(expression)
    matches = object.findall(chromo)
    print(matches[0])

    with open("translate4.txt","r") as translate:
        amino = translate.read()
        amino = amino.split("\n")
        codons = {}
        for codon in amino:
            if codon != "":
                codon = codon.split("\t")
                codons[codon[0]] = codon[1]
        start = 0
        end = 3
        nucleotides = list(mRNA_sequence)
        translation =""
        while end <= len(nucleotides):
            codon = mRNA_sequence[start:end]
            translation = translation + codons.get(codon)
            start = start + 3
            end = end + 3



