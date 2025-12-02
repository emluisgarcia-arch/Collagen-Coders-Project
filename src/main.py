from Bio import SeqIO
import glob
import os

BASE = os.path.dirname(__file__)

ref_file = os.path.join(BASE, "../fasta_data/COL2A1RefSequenceExample.fasta")
disease_folder = os.path.join(BASE, "../fasta_data/")

# Known pathogenic sites for COL2A1
PATHOGENIC_MAP = {
    "Stickler Syndrome": [216, 873, 219, 1300],
    "Achondrogenesis Type 2": [711, 1170, 504],
    "Spondyloepiphyseal Dysplasia": [421, 1504],
    "Osteoarthritis": [1091, 780],
}


def read_sequence(filepath):
    try:
        for record in SeqIO.parse(filepath, "fasta"):
            return str(record.seq)
    except FileNotFoundError:
        print(f"File not found: {filepath}")
    return ""


def compare_sequences(ref_seq, dis_seq):
    differences = []
    length = min(len(ref_seq), len(dis_seq))

    for i in range(length):
        if ref_seq[i] != dis_seq[i]:
            differences.append((i + 1, ref_seq[i], dis_seq[i]))

    if len(ref_seq) != len(dis_seq):
        differences.append(("Length difference", len(ref_seq), len(dis_seq)))

    return differences


def predict_disease(differences):
    matched = []

    # scan all differences
    for pos, ref, var in differences:
        if not isinstance(pos, int):
            continue

        # check if pos is in ANY disease list
        for disease, sites in PATHOGENIC_MAP.items():
            if pos in sites:
                matched.append((disease, pos))

    if not matched:
        return ["Unknown Variant"]

    return matched


# ------------ MAIN ------------

ref_seq = read_sequence(ref_file)

for file in glob.glob(disease_folder + "*.fasta"):

    dis_seq = read_sequence(file)

    differences = compare_sequences(ref_seq, dis_seq)
    predicted = predict_disease(differences)

    filename = os.path.basename(file)

    print("\n==========================")
    print(f"Comparing file: {filename}")
    print("==========================")

    for pos, ref, var in differences:
        print(f"Position {pos}: {ref} -> {var}")

    print("\nPredictions:")
    for d in predicted:
        print("  ", d)

