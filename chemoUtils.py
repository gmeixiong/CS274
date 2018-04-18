# Gerry Meixiong
# CS274 Project 4 chemoUtils.py

import pandas as pd

class chemoUtils:
	# Initialize dictionaries and list for use by tanimoto.py and pvalue.py
	def __init__(self):
		self.fpsList = []
		self.fpsDict = {}
		self.targets = {}
		self.ligands = {}

	# Given a csv file name, returns the matrix values read from pandas
	def matrix(self, file):
		return pd.read_csv(file).values

	# Computes tanimoto from two fingerprints. Fingerprints are sets of keys.
	def computeTanimoto(self, fp1, fp2):
		return len(fp1&fp2)*1.0 / len(fp1|fp2)

	# Add fingerprints to list or dictionary. Each drug has a corresponding element or key.
	# The value is the drug's fingerprint.
	def addFps(self, drugs, listFlag):
		for row in drugs:
			name = row[0]
			fingerprint = set(row[2].split())
			if listFlag:
				self.fpsList.append((name, fingerprint))
			else:
				self.fpsDict[name] = fingerprint
		return

	# Add targets to dictionary. If ligandFlag is true, keys are proteins and values are the drugs/ligands 
	# associated. Otherwise, keys are drugs and values are targets.
	def addTargets(self, targets, ligandFlag):
		for row in targets:
			if ligandFlag:
				key = row[1]
				val = row[0]
				d = self.ligands
			else:
				key = row[0]
				val = row[2]
				d = self.targets
			value_list = d.get(key, [])
			d[key] = value_list + [val]
		return