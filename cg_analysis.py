from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib.pyplot as plt

FILE_NAME = "reads_for_analysis.fastq"

cg_fractions = []

with open(FILE_NAME, "r") as handle:
    for title, seq, qual in FastqGeneralIterator(handle):
        cg_fractions.append(((seq.count("C") + seq.count("G")) / len(seq)) * 100)

print(sorted(cg_fractions))

plt.hist(cg_fractions, bins=20, edgecolor="black")
plt.xlabel("C/G Content (%)")
plt.ylabel("Number of Reads")
plt.title("Distribution of C/G Content in Reads")
plt.show()
