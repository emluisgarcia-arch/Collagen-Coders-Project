from Bio import SeqIO


# Defines the files being compared.
ref_file = "COL2A1RefSequenceExample.fasta"
dis_files = "SticklerExample.fasta"

#Reads the sequence from a FASTA file and returns it as a string.
def read_sequence(filepath):
    try:
        for record in SeqIO.parse(filepath, "fasta"):
            return str(record.seq)
    except FileNotFoundError:
        print(f"File not found at '{filepath}'.")
    return ""

#Compares two sequences and returns a list called differences showing its position in each sequence.
def compare_sequences(ref_seq, var_seq):
    differences = []
    
    ref_len = len(ref_seq)
    var_len = len(var_seq)
    length = min(ref_len, var_len)

    for i in range(length):
        if ref_seq[i] != var_seq[i]:
            differences.append((i + 1, ref_seq[i], var_seq[i]))
    if ref_len != var_len:
        differences.append(("Length difference between sequences", str(ref_len), str(var_len))) 

    return differences

# Predicts a COL2A1-related disease.
def predict_disease(differences):
    for pos, ref_len, var_len in differences:
        if pos == "Length difference between sequences":
            return "Achondrogenesis Type 2 (ACH2)"

    for pos, ref, var in differences:
        if isinstance(pos, int):
            if ref == "G" and var != "G":
                return f"Stickler Syndrome"
            
            #Osteoarthritis and Spondyloepiphyseal need implemntation.
    if differences:
        return "Unknown Variant"

    return "Identical Sequences"
ref_seq = read_sequence(ref_file)
dis_seq = read_sequence(dis_files)
different = compare_sequences(ref_seq, dis_seq)
print(compare_sequences(ref_seq, dis_seq), predict_disease(different))