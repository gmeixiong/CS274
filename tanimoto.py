# Gerry Meixiong
# CS274 Project 4 Tanimoto

import sys
from chemoUtils import chemoUtils

# Parses sys.argv for the proper file names for drugs, targets, and output
def parseArgs(args):
	if len(args) != 4:
		print 'Input format error: python tanimoto.py drugs.csv targets.csv output.csv'
		exit()
	drugs = args[1]
	targets = args[2]
	output = args[3]
	return drugs, targets, output

# Given two lists of targets, return 1 if they share a target, 0 otherwise.
def shareTarget(t1, t2):
	for key in t1:
		if key in t2:
			return 1
	return 0

# Writes each pairwise drug to the output file with the tanimoto score and if they share 
# a target. Only writes lines when both drugs are found in drugs.csv. 
def writeToOutput(cu, fps, targets, output_file):
	f = open(output_file, "w")
	for i in range(0, len(fps)-1):
		for j in range(i+1, len(fps)):
			drug1 = fps[i][0]
			drug2 = fps[j][0]
			if drug1 in targets and drug2 in targets:
				tanimoto = cu.computeTanimoto(fps[i][1], fps[j][1])
				share = shareTarget(targets[drug1], targets[drug2])
				f.write(','.join([drug1, drug2, '{0:6f}'.format(tanimoto), str(share)])+'\n')
	f.close()

# Parses drugs, targets, and output file names. Creates a list of fingerprints from drugs.csv.
# Creates a dictionary of targets from targets.csv. Writes to output.
def main():
	drugs, targets, output = parseArgs(sys.argv)
	cu = chemoUtils()
	cu.addFps(cu.matrix(drugs), True)
	cu.addTargets(cu.matrix(targets), False)
	writeToOutput(cu, cu.fpsList, cu.targets, output)

if __name__ == "__main__":
	main()
