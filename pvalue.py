# Gerry Meixiong
# CS274 Project 4 Pvalue

import argparse
import random
from chemoUtils import chemoUtils

cutoff = 0.5

# Uses argparse to parse command line.
def parseArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', type=int, default=100)
	parser.add_argument('-r', type=int, default=214)
	parser.add_argument('drugs')
	parser.add_argument('targets')
	parser.add_argument('proteinA')
	parser.add_argument('proteinB')
	return parser.parse_args()

# Computes the Tsummary for two sets of ligands. For each pairwise combination of the
# ligand sets, compute the tanimoto score. Add the tanimoto score if it is greater than the cutoff.
def Tsumm(cu, fps, ligandSetA, ligandSetB):
	T = 0
	for ligandA in ligandSetA:
		for ligandB in ligandSetB:
			tanimoto = cu.computeTanimoto(fps[ligandA], fps[ligandB])
			if tanimoto > cutoff:
				T += tanimoto
	return T

# Runs n iterations to compute the bootstrap p-value for the given Tsummary. Calculates tanimoto
# between two random samples of proper length n times. 
def runIterations(Tsum, na, nb, cu, n):
	numGreater = 0
	for i in range(n):
		setA = random.sample(cu.fpsDict, na)
		setB = random.sample(cu.fpsDict, nb)
		randT = Tsumm(cu, cu.fpsDict, setA, setB)
		if randT > Tsum:
			numGreater += 1
	return numGreater*1.0 / n

# Computes bootstrap pvalue of two given proteins. This function is imported by networkgen.py. 
def computeBootstrap(drugs, targets, proteinA, proteinB, r, n):
	random.seed(r)
	cu = chemoUtils()
	cu.addFps(cu.matrix(drugs), False)
	cu.addTargets(cu.matrix(targets), True)
	Tsummary = Tsumm(cu, cu.fpsDict, cu.ligands[proteinA], cu.ligands[proteinB])
	na = len(cu.ligands[proteinA])
	nb = len(cu.ligands[proteinB])
	return runIterations(Tsummary, na, nb, cu, n)

# Main method which parses the args and prints the bootstrap value.
def main():
	args = parseArgs()
	print computeBootstrap(args.drugs, args.targets, args.proteinA, args.proteinB, args.r, args.n)

if __name__ == "__main__":
	main()
