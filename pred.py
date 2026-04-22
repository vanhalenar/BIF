import os

from Bio import Phylo, AlignIO

tree = Phylo.read("tree.tre", "newick")

tree.rooted = True

leaves : dict = {clade: tree.distance(clade) for clade in tree.get_terminals()}

minimum = min(leaves.values())
maximum = max(leaves.values())

leaves = {key: ((value - minimum)/(maximum-minimum)) for key, value in leaves.items()}
leaves = {clade.name: value for clade, value in leaves.items()}
print(leaves)

msa = AlignIO.read("msa.fasta", "fasta")
alignment_length = msa.get_alignment_length()

aminoacids = "ARNDCQEGHILKMFPSTWYV"

conservation_scores = [{x: 0 for x in aminoacids} for _ in range(alignment_length)]

for i in range(alignment_length):
    column = msa[:, i]
    for a in aminoacids:
        score_u = 0
        score_d = 0
        for index in range(len(column)):
            val = column[index]
            value = 1 if val == a else 0
            score_u += value * leaves[msa[index].id]
            score_d += leaves[msa[index].id]
        score = score_u / score_d
        conservation_scores[i][a] = score

print("conservation scores on query: ", conservation_scores[0])
print("conservation scores on number 22: ", conservation_scores[21])

hydrophobicity = {}
polarity = {}
helix_frequency = {}

pairs = ['A/L', 'R/K', 'N/M', 'D/F', 'C/P', 'Q/S', 'E/T', 'G/W', 'H/Y', 'I/V']

properties = {
    "Hydrophobicity": hydrophobicity,
    "Polarity": polarity,
    "Helix frequency": helix_frequency
}

with open("aaindex.txt", "r") as f:
    current = None
    first_row_done = False
    for line in f:
        line = line.strip()
        if line in properties:
            current = properties[line]
            first_row_done = False
        elif line.startswith("A/L"):
            continue
        elif line and current is not None:
            values = line.split()
            if not first_row_done:
                for pair, val in zip(pairs, values):
                    aa1, aa2 = pair.split('/')
                    current[aa1] = float(val)
                first_row_done = True
            else:
                for pair, val in zip(pairs, values):
                    aa1, aa2 = pair.split('/')
                    current[aa2] = float(val)

# min-max normalize each property
for prop in [hydrophobicity, polarity, helix_frequency]:
    mn, mx = min(prop.values()), max(prop.values())
    for aa in prop:
        prop[aa] = (prop[aa] - mn) / (mx - mn)

print("hydrophobicity: ", hydrophobicity)
print("polarity: ", polarity)
print("helix frequency: ", helix_frequency)

query = msa[0]

results = [[0] * alignment_length for _ in range(len(aminoacids))]

for i in range(alignment_length):
    column = msa[:, i]
    original_aminoacid = query.seq[i]
    if original_aminoacid == "-":
        for row in results:
            row[i] = "-"
        continue
    for new_aminoacid in aminoacids:
        original_conservation = conservation_scores[i][original_aminoacid]
        new_conservation = conservation_scores[i][new_aminoacid]
        diff = (original_conservation*hydrophobicity[original_aminoacid] - new_conservation*hydrophobicity[new_aminoacid]) + (original_conservation*polarity[original_aminoacid] - new_conservation*polarity[new_aminoacid]) + (original_conservation*helix_frequency[original_aminoacid] - new_conservation*helix_frequency[new_aminoacid])
        results[aminoacids.index(new_aminoacid)][i] = diff


# output


if os.path.isfile("output.csv"):
    os.remove("output.csv")

with open("output.csv", "w") as f:
    f.write("AA,")
    f.write(",".join(str(i) for i in range(1, alignment_length+1)))
    f.write("\n")

    for i, a in enumerate(aminoacids):
        f.write(a + ",")
        f.write(",".join(str(v) for v in results[i]))
        f.write("\n")

# What's the original amino acid at position 22 (index 21)?                                                                                                                                        
print(f"Query amino acid at position 22: {query.seq[21]}")