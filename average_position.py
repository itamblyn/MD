#! /usr/bin/env python

import sys

print 'talk to Isaac about this...'

import numpy

inputFile = open (sys.argv[1], 'r')
number_of_particles = int(inputFile.readline())
AVERAGE_array = numpy.zeros((number_of_particles,3), dtype=numpy.Float)
lattice_constant = float(sys.argv[2])*0.529177

line_counter = 0

file = inputFile.readlines()
number_of_snapshots = len(file)/(number_of_particles + 2) + 1

for line in file:
     line_split = line.split()
     if line_split[0] == 'H':
          AVERAGE_array[line_counter%number_of_particles][0] += float(line_split[1])
          AVERAGE_array[line_counter%number_of_particles][1] += float(line_split[2])
          AVERAGE_array[line_counter%number_of_particles][2] += float(line_split[3])
          line_counter += 1

inputFile.close()

AVERAGE_array /= float(number_of_snapshots)

outputFile = open('AVERAGE.xyz', 'w')

outputFile.write(str(number_of_particles) + '\n')
outputFile.write('Average of ' + str(number_of_snapshots) + ' snapshots, alat(Bohr) = ' + str(lattice_constant/0.529177) + '\n')


for triple in AVERAGE_array: 
     outputFile.write('Li ')
     for coordinate in triple: 
          pbc_coordinate = coordinate - (int( coordinate/lattice_constant + number_of_snapshots + 0.5) - number_of_snapshots)*lattice_constant
          outputFile.write(str(pbc_coordinate) + ' ')
     outputFile.write('\n')

outputFile.close()
