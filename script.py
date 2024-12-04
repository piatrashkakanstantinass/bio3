from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW, NCBIXML
import csv

FILE_NAME = "reads_for_analysis.fastq"

min_ascii = 999999
max_ascii = -99999
cg_fractions = []


def calc_fraction(seq: str):
    return ((seq.count("C") + seq.count("G")) / len(seq)) * 100


with open(FILE_NAME, "r") as handle:
    for title, seq, qual in FastqGeneralIterator(handle):
        cg_fractions.append(calc_fraction(seq))
        for char in qual:
            ascii_value = ord(char)
            if ascii_value < min_ascii:
                min_ascii = ascii_value
            if ascii_value > max_ascii:
                max_ascii = ascii_value

print(f"Smallest ASCII code: {min_ascii} ('{chr(min_ascii)}')")
print(f"Largest ASCII code: {max_ascii} ('{chr(max_ascii)}')\n")


def guess_encoding(min_ascii, max_ascii):
    possible_encodings = []

    # https://en.wikipedia.org/wiki/FASTQ_format
    encodings = [
        {
            "name": "Sanger Phred+33",
            "min_ascii": 33,
            "max_ascii": 73,
        },
        {
            "name": "Illumina 1.8+ Phred+33",
            "min_ascii": 33,
            "max_ascii": 74,
        },
        {
            "name": "Illumina 1.3+ Phred+64",
            "min_ascii": 64,
            "max_ascii": 104,
        },
        {
            "name": "Illumina 1.5+ Phred+64",
            "min_ascii": 67,
            "max_ascii": 105,
        },
        {
            "name": "Solexa Solexa+64",
            "min_ascii": 59,
            "max_ascii": 104,
        },
    ]

    for enc in encodings:
        if min_ascii >= enc["min_ascii"] and max_ascii <= enc["max_ascii"]:
            possible_encodings.append(enc["name"])

    return possible_encodings


possible_encodings = guess_encoding(min_ascii, max_ascii)

if possible_encodings:
    print("Possible encodings based on ASCII range:")
    for encoding in possible_encodings:
        print(f"- {encoding}")
else:
    print("No standard encoding could be determined from the ASCII range.")

plt.hist(cg_fractions, bins=20, edgecolor="black")
plt.xlabel("C/G Content (%)")
plt.ylabel("Number of Reads")
plt.title("Distribution of C/G Content in Reads")
plt.savefig("graphics.png")

PEAKS = [35, 55]

picked_sequences = []
picked_count = [0, 0]
result_table = []

with open(FILE_NAME, "r") as handle:
    for title, seq, qual in FastqGeneralIterator(handle):
        cg_fraction = calc_fraction(seq)
        for index, peak in enumerate(PEAKS):
            if picked_count[index] >= 5:
                continue
            if abs(cg_fraction - peak) < 0.1:
                picked_sequences.append([title, seq])
                picked_count[index] += 1
                break

for pick in picked_sequences:
    res = NCBIWWW.qblast(
        "blastn", "nt", pick[1], entrez_query="bacteria[organism]", hitlist_size=1
    )
    record = NCBIXML.read(res)
    if record.alignments:
        species_info = record.alignments[0].hit_def
        species = species_info.split("[")[-1].rstrip("]")
    else:
        species = "No match found"
    result_table.append((pick[0], species))

with open("blast_results.csv", "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Read ID", "Species"])
    for read_id, species in result_table:
        csvwriter.writerow([read_id, species])
