from Bio import SeqIO
import argparse

#Reads the sequence from a FASTA file and returns it as a string.
def read_sequence(filepath):
    try:
        for record in SeqIO.parse(filepath, "fasta"):
            return str(record.seq)
    except FileNotFoundError:
        print(f"File not found at '{filepath}'.")
    return ""

#Compares two sequences and returns a list called differences showing its position in each sequence.
def compare_sequences(ref_seq, dis_seq):
    differences = []
    
    ref_len = len(ref_seq)
    var_len = len(dis_seq)
    length = min(ref_len, var_len)

    for i in range(length):
        if ref_seq[i] != dis_seq[i]:
            differences.append((i + 1, ref_seq[i], dis_seq[i]))
    if ref_len != var_len:
        differences.append(("Length difference between sequences", str(ref_len), str(var_len))) 

    return differences

def trimmed_diff(differences):
    if len(differences) > 3:
        return differences[:3] + differences[-1:]
    return differences


# Predicts a COL2A1-related disease.
def predict_disease(differences):
    for pos, ref, dis in differences:
        if isinstance(pos, int):
            if pos == 1000 and ref != dis:
                return "\nStickler Syndrome"
            elif pos == 1500 and ref != dis:
                return "\nOsteoarthritis"
            elif pos == 501 and ref != dis:
                return "\nAchondrogenesis"
            elif pos == 621 and ref != dis:
                return "\nSpondyloepiphyseal"
    if differences:
        return "\nUnknown Variant"

    return "Identical Sequences"

def main():
    parser = argparse.ArgumentParser(
        description="Compares reference sequence to disease sequences to predict COL2A1-related diseases."
    )
    parser.add_argument(
        "-r", "--ref", required=True, help="Reference FASTA file"
    )
    parser.add_argument(
        "-d", "--disease", nargs="+", required=True,
        help="One or more disease FASTA files to compare against the reference"
    )
    args = parser.parse_args()

    ref_seq = read_sequence(args.ref)

    for file in args.disease:
        dis_seq = read_sequence(file)
        print(f'Comparing: {file}')
        compare = compare_sequences(ref_seq, dis_seq)
        differences = trimmed_diff(compare)
        print(differences, predict_disease(compare))


if __name__ == "__main__":
    main()