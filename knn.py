# Gerry Meixiong
# CS274 Project 2 knn
from __future__ import division
import sys
import math
import random
import numpy as np
import pandas as pd

# Can easily swap which label we call positive
pos = 'ALL'
neg = 'AML'

# Given n and a partition number, separates into n random groups. Numbers are
# first shuffled with random.shuffle, then distributed as evenly as possible.
def partition(n, total):
	indexes = [i for i in range(total)]
	random.shuffle(indexes)
	return [[indexes[i] for i in range(j, total, n)] for j in range(n)]

# Given a list of indexes, coalesce the indexed groups from posGroups and
# negGroups. Returns the end grouping of positive and negative indexes.
def getIndexes(indexList, posGroups, negGroups):
	posIndexes, negIndexes = [], []
	for index in indexList:
		posIndexes += posGroups[index]
		negIndexes += negGroups[index]
	return posIndexes, negIndexes

# Given corresponding indexes for the positive and negative files, return
# a list of all desired rows with in order of index, with positive first.
def getSet(posIndexes, negIndexes, pos_file, neg_file):
	dataSet = []
	for index in posIndexes:
		dataSet.append(pos_file[:, index])
	for index in negIndexes:
		dataSet.append(neg_file[:, index])
	return dataSet

# Finds the euclidean distance between the candidate and test points.
# Both points are np.arrays which makes operations much more efficient.
def euclideanDistance(candidatePoint, testPoint):
	diffs = candidatePoint - testPoint
	sqDiffs = diffs**2
	sumOfSqs = np.sum(sqDiffs)
	return np.sqrt(sumOfSqs)

# Returns a list of k integers representing the indexes of the k nearest
# neighbor. Euclidean distance is computed between the target point and 
# each point in the training set. The nearest k are returned. 
def kNearestNeighbors(targetPoint, trainingSet, k):	
	candidatePoints = []
	for index in range(len(trainingSet)):
		distance = euclideanDistance(targetPoint, trainingSet[index])
		candidatePoints.append((index, distance))
	candidatePoints.sort(key=lambda x:x[1])
	neighbors = []
	for index in range(k):
		neighbors.append(candidatePoints[index][0])
	return neighbors

# Given the indexes representing the k neighbors, count the fraction of
# neighbors that are positive. If this is at least p, return a positive
# label. Otherwise, return a negative label.
def getPrediction(neighborIndexes, trainLabels, p, k):
	numPos = 0
	for index in neighborIndexes:
		if trainLabels[index] == pos:
			numPos += 1
	if numPos / k >= p:
		return pos
	else:
		return neg

# Given the test set, returns the predicted labels for each test point based
# on finding the k nearest neighbors and counting the proportion higher than
# p. Look at kNearestNeighbors and getPrediction functions. 
def getPredictions(testSet, trainSet, trainLabels, p, k):
	predLabels = []
	for testPoint in testSet:
		neighborIndexes = kNearestNeighbors(testPoint, trainSet, k)
		predLabels.append(getPrediction(neighborIndexes, trainLabels, p, k))
	return predLabels

# Given the predicted labels and the true labels, return the number of
# true positives, false positives, true negatives, and false negatives.
def getPandN(predLabels, trueLabels):
	TP = FP = TN = FN = 0
	for index in range(len(trueLabels)):
		if trueLabels[index] == pos:
			if predLabels[index] == pos:
				TP += 1
			else:
				FN += 1
		else:
			if predLabels[index] == pos:
				FP += 1
			else:
				TN += 1
	return TP, FP, TN, FN

# Writes all necessary information  with proper format to output file 'knn.out'.
def writeOutput(k, p, n, accuracy, sensitivity, specificity):
	output_file = open('knn.out', 'w')
	output_file.write('k: {}\n'.format(k))
	output_file.write('p: {0:.2f}\n'.format(p))
	output_file.write('n: {}\n'.format(n))
	output_file.write('accuracy: {0:.2f}\n'.format(accuracy))
	output_file.write('sensitivity: {0:.2f}\n'.format(sensitivity))
	output_file.write('specificity: {0:.2f}'.format(specificity))
	output_file.close()

if len(sys.argv) != 6:
	print "Correct input format: python knn.py pos neg k p n"  
else:
	pos_file = pd.read_csv(sys.argv[1], sep = '\t', header = None).values
	neg_file = pd.read_csv(sys.argv[2], sep = '\t', header = None).values

	k = int(sys.argv[3])
	p = float(sys.argv[4])
	n = int(sys.argv[5])

	num_pos = len(pos_file[0])
	num_neg = len(neg_file[0])

	labels = [pos] * num_pos + [neg] * num_neg

	nPosGroups = partition(n, num_pos)
	nNegGroups = partition(n, num_neg)

	totalTP = totalFP = totalTN = totalFN = 0
	for index in range(n):
		allIndexes = [i for i in range(n)]
		trainGroupIndexes = [i for i,j in enumerate(allIndexes) if i != index]
		trainPosIndexes, trainNegIndexes = getIndexes(trainGroupIndexes, nPosGroups, nNegGroups)
		testPosIndexes, testNegIndexes = getIndexes([index], nPosGroups, nNegGroups)
		trainSet = getSet(trainPosIndexes, trainNegIndexes, pos_file, neg_file)
		testSet = getSet(testPosIndexes, testNegIndexes, pos_file, neg_file)
		trainLabels = [pos] * len(trainPosIndexes) + [neg] * len(trainNegIndexes)
		testLabels = [pos] * len(testPosIndexes) + [neg] * len(testNegIndexes)
		predLabels = getPredictions(testSet, trainSet, trainLabels, p, k)
		TP, FP, TN, FN = getPandN(predLabels, testLabels)

		totalTP += TP
		totalFP += FP
		totalTN += TN
		totalFN += FN

	sensitivity = totalTP / (totalTP+totalFN)
	specificity = totalTN / (totalTN+totalFP)
	accuracy = (totalTP+totalTN) / (num_pos+num_neg)
	writeOutput(k, p, n, accuracy, sensitivity, specificity)
