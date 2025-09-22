#Import the needed functions from biopython
from Bio import Entrez, SeqIO

#This is the email NCBI will contact if there are any issues.
Entrez.email = "emluisgarcia@gmail.com"

#refseq id of most common COL2A1
refSeq_ID = "NM_001844.4"

#fetch fasta for COL2A1 gene
handle = Entrez.efetch(db = "nucleotide", id = refSeq_ID, rettype = "fasta", retmode = "text")
seqRecord = SeqIO.read(handle, "fasta")
handle.close()

#print seq record
print(seqRecord.format("fasta"))

#Download on to file
with open("COL2A1RefSequence.fasta", "w") as f:
    SeqIO.write(seqRecord, f, "fasta")