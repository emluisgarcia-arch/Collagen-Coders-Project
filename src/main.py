from Bio import SeqIO
import glob

# Path to the reference sequence
ref_file = "COL2A1RefSequenceExample.fasta"

# Path to folder containing all disease FASTA files
disease_folder = "fasta_data/"     # <-- your folder

# Reads a sequence from FASTA file
def read_sequence(filepath):
    try:
        for record in SeqIO.parse(filepath, "fasta"):
            return str(record.seq)
    except FileNotFoundError:
        print(f"File not found: {filepath}")
    return ""

# Compare reference to variant
def compare_sequences(ref_seq, dis_seq):
    differences = []
    length = min(len(ref_seq), len(dis_seq))

    for i in range(length):
        if ref_seq[i] != dis_seq[i]:
            # store index+1 so it's biological numbering
            differences.append((i + 1, ref_seq[i], dis_seq[i]))

    # If the lengths differ, record that too
    if len(ref_seq) != len(dis_seq):
        differences.append(("Length difference", len(ref_seq), len(dis_seq)))

    return differences

# Very simple disease logic (expand later)
def predict_disease(differences):
    for pos, ref, var in differences:
        if isinstance(pos, int):
            if ref == "G" and var != "G":
                return "Stickler Syndrome"
            # add more conditions for other diseases here

    if differences:
        return "Unknown Variant"

    return "Identical Sequences"

# ---------------- MAIN ----------------

ref_seq = read_sequence(ref_file)

# Loop through every FASTA in folder
for file in glob.glob(disease_folder + "*.fasta"):

    dis_seq = read_sequence(file)
    differences = compare_sequences(ref_seq, dis_seq)
    disease = predict_disease(differences)

    filename = file.split("/")[-1]  # print only file name

    print("\n==========================")
    print(f"Comparing file: {filename}")
    print("==========================")
    print("Differences:", differences)
    print("Predicted Disease:", disease)
