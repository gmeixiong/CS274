# Gerry Meixiong
# CS 274 Project 2 kmeans
from __future__ import division
import sys
import math
import random
import numpy as np
import pandas as pd

# Generates list of k random centroids from data. Use np.amax and np.amin to 
# find centroid bounds. Use random.uniform to pick a point within each range. 
def generateRandomCentroids(k, data):
	dataMaxes = np.amax(data, axis = 0)
	dataMins = np.amin(data, axis = 0)
	centroids = []
	for i in range(k):
		new_centroid = []
		for j in range(len(dataMaxes)):
			new_centroid.append(random.uniform(dataMins[j], dataMaxes[j]))
		centroids.append(new_centroid)
	return centroids

# Finds the euclidean distance between two lists of features, gene and centroid.
# gene and centroid are np.arrays which makes operations more efficient.
def euclideanDistance(gene, centroid):
	diffs = gene - centroid
	sqDiffs = diffs**2
	sumOfSqs = np.sum(sqDiffs)
	return np.sqrt(sumOfSqs)

# Computes distance between gene and each centroid. Returns index of closest centroid.
def findClosest(gene, centroids):
	distances = np.asarray([0] * len(centroids))
	for i in range(len(centroids)):
		distances[i] = euclideanDistance(gene, centroids[i])
	return np.argmin(distances)

# Given a list of centroids and which points belong to them (belongs), move
# each centroid to the center of its points. np.mean can't be called on [].
def moveCentroids(centroids, belongs, geneData):
	for i in range(len(centroids)):
		clusterIndexes = np.where(belongs == i)
		cluster = geneData[clusterIndexes]
		if len(cluster) > 0:
			centroids[i] = np.mean(cluster, axis = 0)
	return centroids

# Writes the cluster placement for each row of features or genes to 'kmeans.out'.
def writeOutput(clusterPlacements):
	output_file = open('kmeans.out', 'w')
	for i in range(len(clusterPlacements)):
		output_file.write('{}\t{}\n'.format(i+1, clusterPlacements[i]+1))
	output_file.close()

if len(sys.argv) == 5:
	centroids = pd.read_csv(sys.argv[4], sep = '\t', header = None).values
else:
	centroids = None

if len(sys.argv) < 4 or len(sys.argv) > 5:
	print "Correct input format: python kmeans.py k express.dat max.it [centroids.txt]"
else:
	data = pd.read_csv(sys.argv[2], sep = '\t', header = None).values
	k = int(sys.argv[1])
	maxIters = int(sys.argv[3])
	if centroids is None:
		centroids = np.asarray(generateRandomCentroids(k, data))
	else:
		centroids = centroids[0:k,]

	clusterPlacements = np.asarray([None] * len(data))
	for numIters in range(maxIters+1):
		sameClusters = True
		for i in range(len(data)):
			geneData = data[i]
			closestCluster = findClosest(geneData, centroids)
			if clusterPlacements[i] != closestCluster:
				clusterPlacements[i] = closestCluster
				sameClusters = False
		if sameClusters:
			break
		else:
			centroids = moveCentroids(centroids, clusterPlacements, data)
	writeOutput(clusterPlacements)
	print 'iterations: {}'.format(numIters)

