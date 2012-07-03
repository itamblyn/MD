#! /usr/bin/env python

import sys
import commands
import numpy

print 'usage: ' + sys.argv[0] + ' number_of_bins lattice_constant(bohr) xyz_file'

def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
          if (input_value > 0): i+= 1
          if (input_value < 0): i-= 1
     return i

number_of_bins, lattice_constant, input_filename, snapshot_stride = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], 1 # should be a user input

if snapshot_stride < 1:
    snapshot_stride = 1

command_line_counter = commands.getoutput('wc -l ' + input_filename).split()

if len(command_line_counter) != 2:
    print 'Error determining file size'
else:
    number_of_lines = int(command_line_counter[0])

size_of_box_x, size_of_box_y, size_of_box_z = lattice_constant, lattice_constant, lattice_constant

xyzFile = open(input_filename, 'r')

number_of_particles = int(xyzFile.readline())
number_of_snapshots = number_of_lines/(number_of_particles + 2)         # the 2 is due to xyz format
number_of_neighbours = int(number_of_particles - 1)


###### histgram stuff

min_distance = 0.0
max_distance = lattice_constant/2.0

histogram = numpy.zeros((int(number_of_neighbours), int(number_of_bins)), dtype=numpy.int)
bin_size = (max_distance - min_distance)/float(number_of_bins)


SNAPSHOT_array = numpy.zeros((number_of_particles, 3), dtype=numpy.float)

xyzFile.seek(0)			# return to beginning of file

s = 0
while s < number_of_snapshots:

    distance_array = []

    xyzFile.readline()          # skip number of particles
    xyzFile.readline()          # skip MD step or cell parameters

    for p in range(number_of_particles):
        line = xyzFile.readline()
        SNAPSHOT_array[p][0] = float(line.split()[1])/0.529177
        SNAPSHOT_array[p][1] = float(line.split()[2])/0.529177
        SNAPSHOT_array[p][2] = float(line.split()[3])/0.529177

    # At this point, SNAPSHOT_array is good to go!!!

    for p in range(number_of_particles):        # counts over particles 

         minimum_distances = numpy.zeros(number_of_particles, dtype=numpy.float)

         for o in range(number_of_particles):

              dx = SNAPSHOT_array[p][0] - SNAPSHOT_array[o][0]
              dy = SNAPSHOT_array[p][1] - SNAPSHOT_array[o][1]
              dz = SNAPSHOT_array[p][2] - SNAPSHOT_array[o][2]

              dx -= size_of_box_x*pbc_round(dx/size_of_box_x)
              dy -= size_of_box_y*pbc_round(dy/size_of_box_y)
              dz -= size_of_box_z*pbc_round(dz/size_of_box_z)

              distance = (dx**2 + dy**2 + dz**2)**(0.5)

              minimum_distances[o] = distance

         minimum_distances = numpy.sort(minimum_distances)

         distance_array.append([])
         for i in range(number_of_particles - 1):
              distance_array[-1].append(minimum_distances[i + 1])

    # At this point, distance_array is ready to go

    for neighbour_list in distance_array:
        for neighbour_index in range(number_of_neighbours):
            if neighbour_list[neighbour_index] < max_distance:
                bin = int((neighbour_list[neighbour_index] - min_distance)/bin_size)
                histogram[neighbour_index][bin] += 1

    s += snapshot_stride


outputFile_nn_average = open ('nn_average.hist', 'w')
outputFile_nn_average.write('# bin (Angst), nn1, nn2, ... \n')


for bin_index in range(number_of_bins):

     outputFile_nn_average.write(repr( (bin_index*bin_size + bin_size/2.0 + min_distance)*0.529177 ) + '  ')
     for neighbour_index in range(number_of_neighbours):
          outputFile_nn_average.write(str(histogram[neighbour_index][bin_index]/float(number_of_particles*number_of_snapshots)) + ' ')
     outputFile_nn_average.write('\n')

outputFile_nn_average.close()
