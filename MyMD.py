# Gerry Meixiong
# CS274 Project 3

from __future__ import division
import sys
import math
import numpy as np

# Error with too few or too many arguments
if len(sys.argv) < 3 or len(sys.argv) > 17:
	print "Input format error: python MyMD.py --if input_file --kb [] ..."
	exit()

# Initialize parameters of interest.
input_file = None
kb = 40000.0
kn = 400.0
nbCutoff = 0.50
m = 12.0
dt = 0.001
n = 1000
out = None

# Initialize lists of interest. Every atom will have an element in these lists.
positions = []
velocities = []
b_connections = []
nb_connections = []
b_refs = []
nb_refs = []

# Parses command line input and replaces default values with given parameters.
#
# Input: list of command line arguments following python MyMD.py
# Output: parameters are replaced if given in command line. otherwise set to default.
def parseInput(args):
	global input_file, kb, kn, nbCutoff, m, dt, n, out
	for i in range(0, len(args), 2):
		if args[i].lower() == '--if':
			input_file = args[i+1]
			if not out:
				out = input_file[:-4]
		elif args[i].lower() == '--kb': 
			kb = float(args[i+1])
		elif args[i].lower() == '--kn':
			kn = float(args[i+1])
		elif args[i].lower() == '--nbcutoff':
			nbCutoff = float(args[i+1])
		elif args[i].lower() == '--m':
			m = float(args[i+1])
		elif args[i].lower() == '--dt':
			dt = float(args[i+1])
		elif args[i].lower() == '--n':
			n = int(args[i+1])
		elif args[i].lower() == '--out':
			out = args[i+1]
		else:
			print "Invalid parameter name."
			exit()

# Parses input file and fills the position, velocity, and bond connection lists.
#
# Input: input file, lists for position, velocity, and bond connections. 
# Output: header of the input file to be printed in output rvc file. Each atom
# 		  now has a corresponding element in each list enumerating its initial
#		  position, velocity, as well as the atoms it directly bonds to. 
def parseInputFile(file, positions, velocities, b_connections):
	f = open(file, 'r')
	header = f.readline()
	for line in f:
		atom_info = line.split()
		position = np.asarray(map(float, atom_info[1:4]))
		velocity = np.asarray(map(float, atom_info[4:7]))
		b_connection = map(lambda x: int(x)-1, atom_info[7:])
		positions.append(position)
		velocities.append(velocity)
		b_connections.append(b_connection)
	f.close()
	return header

# Calculates euclidean distance between two positions. Uses numpy.
#
# Input: two numpy arrays which represent two positions.
# Output: the euclidean distance between the positions.
def euclideanDistance(pos1, pos2):
	diffs = pos1 - pos2
	squares = diffs**2
	return np.sqrt(np.sum(squares))

# Finds the non-bonded interacting atoms and reference lengths for every atom.
# Non-bonded interacting atoms are atoms within a certain distance (nbCutoff) that 
# are not already bonded to the atom. 
#
# Input: atom positions, atom bond connections, nbCutoff, empty list 
# 		 of non-bonded connections and reference lengths to be filled.
# Output: Filled lists of non-bonded connections and corresponding reference lengths.
#		  Every atom has an entry in the list. 
def findNbConnections(positions, b_cons, nbCutoff, nb_connections, nb_reference_lengths):
	for i in range(len(positions)):
		nb_neighbors = []
		nb_ref_length = []
		for j in range(len(positions)):
			if i != j and j not in b_cons[i]:
				dist = euclideanDistance(positions[i], positions[j])
				if dist < nbCutoff:
					nb_neighbors.append(j)
					nb_ref_length.append(dist)
		nb_connections.append(nb_neighbors)
		nb_reference_lengths.append(nb_ref_length)

# Calculates reference lengths between each pair of bonded atoms. Fills the given list.
#
# Input: atom positions, atom bonded connections, empty list of reference lengths.
# Output: Filled list of bonded-connection reference lengths. Every atom has an entry in
#		  this list.
def calculateReferences(positions, connections, b_reference_lengths):
	for i in range(len(connections)):
		b_ref_length = []
		for neighbor_index in connections[i]:
			dist = euclideanDistance(positions[i], positions[neighbor_index])
			b_ref_length.append(dist)
		b_reference_lengths.append(b_ref_length)

# Resets the vectors for force, bond potential energy, and non-bond potential energy.
#
# Input: forces list, bond pe list, nonbond pe list. Each element of forces is a numpy
#		 array to allow for easy updating of components. Potential energy lists are both
#		 numpy arrays to allow for summing to find total potential energy.
# Output: zeroed out lists for force, bond and nonbond potential energy.
def resetForcesAndEnergies(forces, bond_p_energies, nbond_p_energies):
	for force in forces:
		force.fill(0)
	bond_p_energies.fill(0)
	nbond_p_energies.fill(0)

# Updates the energy and force lists for a given list of connections and reference angles. 
# This will be called once for bonded and once for nonbonded.
#
# Input: atom positions, atom connections (either bonded or nonbonded), corresponding reference
#		 lengths, zeroed forces list, zeroed potential energy list, corresponding k value.
# Output: Filled forces and potential energy lists. Adds potential energy on an atom from each of 
#		  its connected atoms. Calculates a force magnitude and adds proper force components for
#		  each atom using spring formulas and derivative of potential energy.
def updateEnergyAndForce(positions, connections, refs, forces, potential_energies, k):
	for i in range(len(connections)):
		for j in range(len(connections[i])):
			neighbor = connections[i][j]
			b = euclideanDistance(positions[i], positions[neighbor])
			potential_energies[i] += 0.5 * k * (b-refs[i][j])**2
			force_magnitude = k * (b - refs[i][j])
			forces[i] += force_magnitude * (positions[neighbor] - positions[i]) / b

# Computes kinetic energy given a mass and list of velocities. Uses numpy.
#
# Input: mass value and numpy array of velocities.
# Output: the total kinetic energy of the system. 
def computeKineticEnergy(m, v):
	return np.sum(0.5*m*v**2)

# Writes the proper energy values and time step to the output erg file.
#
# Input: the opened output erg file, timestep, KE, bPE, nbPE, and totalE to write.
# Output: One tab-delimited line with timestep and energies has been written to file.
def writeToErg(file, ts, KE, bPE, nbPE, totalE):
	items = [ts, KE, bPE, nbPE, totalE]
	rounded = [np.around(_, 1) for _ in items]
	stringItems = map(str, rounded)
	writeString = '\t'.join(stringItems)
	file.write(writeString + '\n')

# Writes the rvc values at a given timestep to the rvc output file.
#
# Input: the opened rvc file, positions, velocities, connections, timestep, and energy.
# Output: A line printing the time step and total energy at that level is written. Then
#		  the position, velocity, and bonded atoms of every atom is written on its own line.
def writeToRvc(file, pos, vel, cons, ts, energy):
	if ts:
		timeString = '#At time step ' + str(ts)
		energyString = ',energy = ' + '{0:.3f}'.format(energy) + 'kJ'
		file.write(timeString + energyString + '\n')

	for atom in range(len(positions)):
		position = pos[atom].tolist()
		velocity = vel[atom].tolist()
		bonds = [i+1 for i in cons[atom]]
		rounded = map('{0:.4f}'.format, position + velocity)
		items = [str(atom+1)] + rounded + map(str,bonds)
		writeString = '\t'.join(items)
		if len(bonds) == 0:
			writeString += '\t'
		file.write(writeString + '\n')

# Parse input.
parseInput(sys.argv[1:])
header = parseInputFile(input_file, positions, velocities, b_connections)

# Find all non-bonded connections and calculate reference lengths for all interacting atoms.
findNbConnections(positions, b_connections, nbCutoff, nb_connections, nb_refs)
calculateReferences(positions, b_connections, b_refs)

# Set up velocity, force, and potential energy vectors.
velocities = np.asarray(velocities)
forces = [np.asarray([0.0,0.0,0.0]) for i in range(len(positions))]
b_PEs = np.asarray([0.0] * len(positions))
nb_PEs = np.asarray([0.0] * len(positions))

# Open files for writing. Write initial lines to files before any iterations.
erg_file, rvc_file = out + '_out.erg', out + '_out.rvc'
erg = open(erg_file, 'w')
rvc = open(rvc_file, 'w')
erg.write('# step\tE_k\tE_b\tE_nB\tE_tot\n')
rvc.write(header)
writeToRvc(rvc, positions, velocities, b_connections, 0, 0)
initialEnergy = None

# Begin simulation of n time-steps. 
for time_step in range(1, n+1):

	# Update velocity at dt/2 timestep. Update position at dt timestep.
	for atom_index in range(len(velocities)):
		a = forces[atom_index] / m
		velocities[atom_index] += 0.5 * a * dt
		positions[atom_index] += velocities[atom_index] * dt

	# Reset force and potential energy. Update new force and energy from new positions.
	resetForcesAndEnergies(forces, b_PEs, nb_PEs)
	updateEnergyAndForce(positions, b_connections, b_refs, forces, b_PEs, kb)
	updateEnergyAndForce(positions, nb_connections, nb_refs, forces, nb_PEs, kn)

	# Update velocity at dt timestep. 
	for atom_index in range(len(velocities)):
		a = forces[atom_index] / m
		velocities[atom_index] += 0.5 * a * dt

	# Calculate kinetic energy, bond and nonbond potential energy, and total
	KE = computeKineticEnergy(m, velocities)
	bondPE = np.sum(b_PEs) / 2
	nonbondPE = np.sum(nb_PEs) / 2
	totalEnergy = bondPE + nonbondPE + KE

	# Initialize initialEnergy on first iteration. If energy growth exceeds factor of 10, 
	# molecule is unstable. Print timestep and exit. 
	if not initialEnergy:
		initialEnergy = totalEnergy
	if totalEnergy / initialEnergy >= 10:
		print 'The molecule is unstable. ' + str(totalEnergy) \
			+ ' exceeds the limit. Timestep: ' + str(time_step)
		break

	# Write to erg and rvc files every 10 time steps.
	if time_step % 10 == 0:
		writeToErg(erg, time_step, KE, bondPE, nonbondPE, totalEnergy)
		writeToRvc(rvc, positions, velocities, b_connections, time_step, totalEnergy)

# Close file descriptors opened for writing.
erg.close()
rvc.close()