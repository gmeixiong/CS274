#!/usr/bin/env python
# Gerry Meixiong
# CS 274 Project 1 Code
import sys

# Function to create matchMatrix with NX rows and NY cols by reading from file.
#
# Input: input_file to be read
# Output: NX x NY match matrix
def createMatchMatrix(input_file):
	matchMatrix = []
	for row in range(NX):
		newRow = []
		for col in range(NY):
			r, c, a, b, s = input_file.readline().strip().split()
			newRow.append(float(s))
		matchMatrix.append(newRow)
	return matchMatrix

# Function that initializes the three score matrices and three corresonding
# backtrace matrices. Score entries initialized to 0. Backtrace entries are 
# initialized to empty lists that will store tuples (arrows).
#
# Input: lengths of sequences X and Y
# Output: Initialized score and backtrace matrices
def initScoreMatrices(lenX, lenY):
	M = [[0.0 for col in range(lenY+1)] for row in range(lenX + 1)]
	IX = [[0.0 for col in range(lenY+1)] for row in range(lenX + 1)]
	IY = [[0.0 for col in range(lenY+1)] for row in range(lenX + 1)]
	Mbt = [[[] for col in range(lenY+1)] for row in range(lenX + 1)]
	IXbt = [[[] for col in range(lenY+1)] for row in range(lenX + 1)]
	IYbt = [[[] for col in range(lenY+1)] for row in range(lenX + 1)]
	return M, IX, IY, Mbt, IXbt, IYbt

# Function that draws the correct arrows for a given cell in the M matrix.
# If score is equal to M[i][j], this signals the start of a local alignment
# so no arrow is necessary. M[i][j] should already be filled at this point.
#
# Input: score for the letters given by i and j; our three candidates from
#        equation 2.16; i, and j which indicate the given cell, local.
# Output: n/a. Appends proper tuple (arrow) to cell in  M's backtrace matrix.
def pickDirectionM(score, choices, i, j, local):
	delta = 0.001
	for index in range(len(choices)):
		if abs(choices[index] - M[i][j]) < delta:
			if local and score == M[i][j]:
				continue
			if index == 0:
				Mbt[i][j].append(('Mbt', i-1, j-1))
			elif index == 1:
				Mbt[i][j].append(('IXbt', i-1, j-1))
			elif index == 2:
				Mbt[i][j].append(('IYbt', i-1, j-1))

# Function that draws the correct arrows for a given cell in the IX matrix.
#
# Input: two IX candidates from equation 2.16; i, and j.
# Output: n/a. Appends proper tuple (arrow) to cell in  IX's backtrace matrix.
def pickDirectionIX(choices, i, j):
	delta = 0.001
	for index in range(len(choices)):
		if abs(choices[index] - IX[i][j]) < delta:
			if index == 0:
				IXbt[i][j].append(('Mbt', i-1, j))
			else:
				IXbt[i][j].append(('IXbt', i-1, j))

# Function that draws the correct arrows for a given cell in the IY matrix.
#
# Input: two IY candidates from equation 2.16; i, and j.
# Output: n/a. Appends proper tuple (arrow) to cell in  IY's backtrace matrix.
def pickDirectionIY(choices, i, j):
	delta = 0.001
	for index in range(len(choices)):
		if abs(choices[index] - IY[i][j]) < delta:
			if index == 0:
				IYbt[i][j].append(('Mbt', i, j-1))
			else:
				IYbt[i][j].append(('IYbt', i, j-1))

# Function that populates the 3 scoring matrices from left-right, top-bottom.
# Starts from (1,1) as all the 0-indexed elements are initialized to 0. Uses 
# equations from 2.16 to generate lists of candidates for each scoring matrix.
# Finds the max of these candidates and populates the respective backtrace 
# matrices with the proper tuples (arrows).
#
# Input: n/a. Function is called so the necessary matrices and alphabet
# 		 information are in scope.
# Output: n/a. At the end of the function call, all score matrices and 
#         backtrace matrices are fully populated. 
def populateMatrices():
	for i in range(1, len(seqX) + 1):
		matchIndexI = alphabetX.find(seqX[i-1])
		for j in range(1, len(seqY) + 1):
			matchIndexJ = alphabetY.find(seqY[j-1])
			scoreIJ = matchMatrix[matchIndexI][matchIndexJ]	
			Mchoices = [M[i-1][j-1] + scoreIJ, IX[i-1][j-1] + scoreIJ, IY[i-1][j-1] + scoreIJ]
			IXchoices = [M[i-1][j] - dy, IX[i-1][j] - ey]
			IYchoices = [M[i][j-1] - dx, IY[i][j-1] - ex]
			if local:
				Mchoices.append(0)
				IXchoices.append(0)
				IYchoices.append(0)
			M[i][j] = max(Mchoices)
			IX[i][j] = max(IXchoices)
			IY[i][j] = max(IYchoices)
			pickDirectionM(scoreIJ, Mchoices, i, j, local)
			pickDirectionIX(IXchoices, i, j)
			pickDirectionIY(IYchoices, i, j)

# Function that finds the max score from all of the scoring matrices.
#
# Input: M, IX, IY scoring matrices; boolean indicating if local.
# Output: max score (for either global or local alignment).
def findMaxAlignmentScore(M, IX, IY, local):
	maxScore = 0
	for row in range(1, len(M)):
		for col in range(1, len(M[0])):
			if not local and row < len(M)-1 and col < len(M[0])-1:
				continue
			maxScore = max(M[row][col], IX[row][col], IY[row][col], maxScore)
	return maxScore

# Function to find the indices and respective table that contains maxScore.
#
# Input: maxScore; M, IX, IY scoring matrices; boolean indicating if local.
# Output: Lists of tables, i-indexes, and j-indexes that enumerate all elements
#         which are equivalent to maxScore. I chose to use a list because local
#         alignment may have multiple cells which are equivalent to maxScore. 
def findMaxIndices(maxScore, M, IX, IY, local):
	delta = 0.001
	tables, i, j = [], [], []
	for row in range(1, len(M)):
		for col in range(1, len(M[0])):
			if not local and row < len(M)-1 and col < len(M[0])-1:
				continue
			if abs(M[row][col] - maxScore) < delta:
				tables.append('Mbt')
				i.append(row)
				j.append(col)
			elif abs(IX[row][col] - maxScore) < delta:
				tables.append('IXbt')
				i.append(row)
				j.append(col)
			elif abs(IY[row][col] - maxScore) < delta:
				tables.append('IYbt')
				i.append(row)
				j.append(col)
	return tables, i, j

# Function to backtrace from our max scoring element and generate the proper
# alignment strings. Recursively builds the alignment string by appending 
# characters one at a time. Handles forking of multiple arrows by iterating
# through tuples. Base case where either the row or col index is 0.
#
# Input: input string X and Y; list of alignments to print; table, row, col 
#        indicating which cell to backtrace from.
# Output: n/a. Populated list of max scoring alignments to be printed.
def generateBacktrace(stringX, stringY, alignments, table, row, col):
	if row == 0 or col == 0: # base case
		if (stringX, stringY) not in alignments:
			alignments.append((stringX, stringY))
		return

	if table == 'Mbt':
		stringX = seqX[row-1] + stringX
		stringY = seqY[col-1] + stringY
		if len(Mbt[row][col]) == 0: # no more to backtrace, start of local
			alignments.append((stringX, stringY))
			return
		for nextTable, i, j in Mbt[row][col]:
			if local and (i == 0 or j == 0):
				generateBacktrace(stringX, stringY, alignments, 'Mbt', i, j)
				break
			else:
				generateBacktrace(stringX, stringY, alignments, nextTable, i, j)
	elif table == 'IXbt':
		stringX = seqX[row-1] + stringX
		stringY = '_' + stringY
		for nextTable, i, j in IXbt[row][col]:
			if local and (i == 0 or j == 0):
				generateBacktrace(stringX, stringY, alignments, 'Mbt', i, j)
				break
			else:
				generateBacktrace(stringX, stringY, alignments, nextTable, i, j)
	else:
		stringX = '_' + stringX
		stringY = seqY[col-1] + stringY
		for nextTable, i, j in IYbt[row][col]:
			if local and (i == 0 or j == 0):
				generateBacktrace(stringX, stringY, alignments, 'Mbt', i, j)
				break
			else:
				generateBacktrace(stringX, stringY, alignments, nextTable, i, j)

# Function to remove the top alignments that have leading gaps.
#
# Input: list of top alignments
# Output: list of top alignments with leading gaps removed.
def truncate(topAlignments):
	return [align for align in topAlignments if align[0][0] != '_' and align[1][0] != '_']

# Main function driving the program. Handles input format, file reading, and calls 
# all helper functions that implement the algorithm. Prints winning alignments to 
# the given output file. 
if len(sys.argv) != 3:
	print "Correct input format: python align.py <input file> <output file>"  
else:
	input_file = open(sys.argv[1], 'r')
	seqX = input_file.readline().strip()
	seqY = input_file.readline().strip()
	local = int(input_file.readline())
	dx, ex, dy, ey = map(float, input_file.readline().split())
	NX = int(input_file.readline())
	alphabetX = input_file.readline().strip()
	NY  = int(input_file.readline())
	alphabetY = input_file.readline().strip()
	matchMatrix = createMatchMatrix(input_file)
	input_file.close()

	M, IX, IY, Mbt, IXbt, IYbt = initScoreMatrices(len(seqX), len(seqY))
	populateMatrices()
	maxScore = findMaxAlignmentScore(M, IX, IY, local)
	tables, maxIs, maxJs = findMaxIndices(maxScore, M, IX, IY, local)

	alignments = []
	for table, maxI, maxJ in zip(tables, maxIs, maxJs):
		generateBacktrace('', '', alignments, table, maxI, maxJ)
	alignments = truncate(alignments)
	
	output_file = open(sys.argv[2], 'w')
	output_file.write("{0:.1f}".format(maxScore))
	output_file.write('\n')
	for alignmentX, alignmentY in alignments:
		output_file.write('\n')
		output_file.write(alignmentX)
		output_file.write('\n')
		output_file.write(alignmentY)
		output_file.write('\n')
	output_file.close()