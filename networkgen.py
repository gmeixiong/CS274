# Gerry Meixiong
# CS274 Project 4 Networkgen

import argparse
from pvalue import computeBootstrap

# Parameters to use for computing bootstrap.
bootstrap_iterations = 100
random_seed = 214

# Uses argparse to parse csv file names for drugs, targets, and protein_nodes.
def parseArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument('drugs')
	parser.add_argument('targets')
	parser.add_argument('proteins')
	return parser.parse_args()

# Writes to .nodeAttr files given a header, protein names, and labels.
# Labels will be indications or ids. 
def writeToFile(outputFile, header, proteins, labels):
	output = open(outputFile, 'w')
	output.write(header + '\n')
	for i in range(len(proteins)):
		output.write(proteins[i] + ' = ' + labels[i] + '\n')
	output.close()

# Parses the protein_nodes csv file for the protein id, protein name, and indications.
# Appends this information the the lists ids, names, and indications.
def parseInput(protein_csv, ids, names, indications):
	protein_file = open(protein_csv, 'r')
	protein_file.readline()
	for row in protein_file:
		info = row.strip().split(',')
		node_id = info[0]
		node_name = info[1]
		indication = info[2]
		ids.append(node_id)
		names.append(node_name)
		indications.append(indicationi)
	protein_file.close()
	return ids, names, indications

# Generates .sif file by writing the pairwise proteins in protein_nodes.csv and computes
# the bootstrap p value between them. Writes the edge to the sif file if p value < 0.05.
def generateSIF(ids, drug_file, target_file):
	output = open('network.sif', 'w')
	for i in range(0, len(ids)-1):
		for j in range(i+1, len(ids)):
			p_bootstrap = computeBootstrap(drug_file, target_file, \
				ids[i], ids[j], random_seed, bootstrap_iterations)
			if p_bootstrap <= 0.05:
				output.write(ids[i] + ' edge ' + ids[j] + '\n')
	output.close()

# Main method which parses the command line, input, and generates the appropriate .nodeAttr
# and .sif files. 
def main():
	args = parseArgs()
	ids, names, indications = parseInput(args.proteins, [], [], [])
	writeToFile('name.nodeAttr', 'name', ids, names)
	writeToFile('indication.nodeAttr', 'indication', ids, indications)
	generateSIF(sorted(ids), args.drugs, args.targets)

if __name__ == "__main__":
	main()
