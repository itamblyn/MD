#! /usr/bin/env python

import sys
import commands
import numpy

print 'usage: ' + sys.argv[0] + ' number_of_bins lx ly lz (bohr) xyz_file'

def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
          if (input_value > 0): i+= 1
          if (input_value < 0): i-= 1
     return i

number_of_bins = int(sys.argv[1])

size_of_box_x = float(sys.argv[2])
size_of_box_y = float(sys.argv[3])
size_of_box_z = float(sys.argv[4])

input_filename = sys.argv[5]

snapshot_stride = 1 # should be a user input

if snapshot_stride < 1:
    snapshot_stride = 1

command_line_counter = commands.getoutput('wc -l ' + input_filename).split()

if len(command_line_counter) != 2:
    print 'Error determining file size'
else:
    number_of_lines = int(command_line_counter[0])

#size_of_box_x, size_of_box_y, size_of_box_z = lattice_constant, lattice_constant, lattice_constant

xyzFile = open(input_filename, 'r')

number_of_particles = int(xyzFile.readline())
number_of_snapshots = number_of_lines/(number_of_particles + 2)         # the 2 is due to xyz format
number_of_neighbours = int(number_of_particles - 1)

outputFile = open('TRAJEC.cnn', 'w')

outputFile.write('# comment =  \n')
outputFile.write('# a = ' + str(size_of_box_x) + '\n')
outputFile.write('# b = ' + str(size_of_box_y) + '\n')
outputFile.write('# c = ' + str(size_of_box_z) + '\n')
outputFile.write('# number_of_particles = '  + str(number_of_particles) + '\n')
outputFile.write('# number_of_neighbours = ' + str(number_of_neighbours) + '\n')
outputFile.write('#\n')
outputFile.write('#\n')
outputFile.write('#\n')
outputFile.write('# units = bohr\n')

###### histgram stuff

min_distance = 0.0
max_distance = max(size_of_box_x, size_of_box_y, size_of_box_z)/2.0

histogram = numpy.zeros((int(number_of_neighbours), int(number_of_bins)), dtype=numpy.int)
bin_size = (max_distance - min_distance)/float(number_of_bins)


SNAPSHOT_array = numpy.zeros((number_of_particles, 3), dtype=numpy.float)

xyzFile.seek(0)			# return to beginning of file

for s in range(number_of_snapshots):

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
         
         outputFile.write('Li' + ' % .8e' % SNAPSHOT_array[p][0] + ' % .8e' % SNAPSHOT_array[p][1] + ' % .8e' % SNAPSHOT_array[p][2])

         neighbour_distances = [] 

         for o in range(number_of_particles):

              dx = SNAPSHOT_array[p][0] - SNAPSHOT_array[o][0]
              dy = SNAPSHOT_array[p][1] - SNAPSHOT_array[o][1]
              dz = SNAPSHOT_array[p][2] - SNAPSHOT_array[o][2]

              dx -= size_of_box_x*pbc_round(dx/size_of_box_x)
              dy -= size_of_box_y*pbc_round(dy/size_of_box_y)
              dz -= size_of_box_z*pbc_round(dz/size_of_box_z)

              distance = (dx**2 + dy**2 + dz**2)**(0.5)

              neighbour_distances.append([distance, o])

         list.sort(neighbour_distances)

#         if s > 0:

#             for monkey in neighbour_distances:
#                 print monkey

         distance_array.append([])

         for i in range(number_of_neighbours):
              distance_array[-1].append(neighbour_distances[i + 1][0])
              outputFile.write(' ' + str.rjust(str(neighbour_distances[i + 1][1]), 3) ) 

         outputFile.write('\n')

    # At this point, distance_array is ready to go

    for distance_list in distance_array:
        for neighbour_index in range(number_of_neighbours):
            if distance_list[neighbour_index] < max_distance:
                bin = int((distance_list[neighbour_index] - min_distance)/bin_size)
                histogram[neighbour_index][bin] += 1

outputFile_nn_average = open ('nn_average.hist', 'w')
outputFile_nn_average.write('# bin (Angst), nn1, nn2, ... \n')


for bin_index in range(number_of_bins):

     outputFile_nn_average.write(repr( (bin_index*bin_size + bin_size/2.0 + min_distance)*0.529177 ) + '  ')
     for neighbour_index in range(number_of_neighbours):
          outputFile_nn_average.write(str(histogram[neighbour_index][bin_index]/float(number_of_particles*number_of_snapshots)) + ' ')
     outputFile_nn_average.write('\n')

outputFile_nn_average.close()
