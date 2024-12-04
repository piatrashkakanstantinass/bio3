from Bio.SeqIO.QualityIO import FastqGeneralIterator

FILE_NAME = "reads_for_analysis.fastq"

min_ascii = 999999
max_ascii = -99999

with open(FILE_NAME, "r") as handle:
    for title, seq, qual in FastqGeneralIterator(handle):
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
