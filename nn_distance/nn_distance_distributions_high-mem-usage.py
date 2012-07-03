#! /usr/bin/env python

import sys

print 'usage: ' + sys.argv[0] + ' number_of_bins lattice_constant(bohr) xyz_file'

import numpy

def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
          if (input_value > 0): i+= 1
          if (input_value < 0): i-= 1
     return i

number_of_bins = int(sys.argv[1])
lattice_constant = float(sys.argv[2])

size_of_box_x = lattice_constant
size_of_box_y = lattice_constant
size_of_box_z = lattice_constant

xyzFile = open (sys.argv[3],'r')

number_of_particles = int(xyzFile.readline())

TRAJECTORY_array = []

for line in xyzFile:
    if len(line.split()) == 4:

         position_vector = [float(line.split()[1])/0.529177, float(line.split()[2])/0.529177,float(line.split()[3])/0.529177]
         TRAJECTORY_array.append(position_vector)

xyzFile.close()

number_of_snapshots = len(TRAJECTORY_array)/float(number_of_particles)

global_min = 1E6
global_max = 0.0

distance_array = []

s = 0 # counts over snapshots

while s < number_of_snapshots:

     p = 0 # counts over particles
     
     while p < number_of_particles:
     
          o = 0 # counts over other particles
          
          minimum_distances = numpy.zeros(number_of_particles, dtype=numpy.float)
          
          while o < number_of_particles:

               dx = TRAJECTORY_array[s*number_of_particles + p][0] - TRAJECTORY_array[s*number_of_particles + o][0] 
               dy = TRAJECTORY_array[s*number_of_particles + p][1] - TRAJECTORY_array[s*number_of_particles + o][1]
               dz = TRAJECTORY_array[s*number_of_particles + p][2] - TRAJECTORY_array[s*number_of_particles + o][2]
               
               dx -= size_of_box_x*pbc_round(dx/size_of_box_x)
               dy -= size_of_box_y*pbc_round(dy/size_of_box_y)
               dz -= size_of_box_z*pbc_round(dz/size_of_box_z)
               
               distance = (dx**2 + dy**2 + dz**2)**(0.5)

               minimum_distances[o] = distance
               
               o +=1
          
	  sorted = numpy.sort(minimum_distances)

          if (sorted[1] < global_min): global_min = sorted[1]
          if (sorted[number_of_particles - 1] > global_max): global_max = sorted[number_of_particles - 1]
          distance_array.append([])
          for i in range(number_of_particles - 1):
               distance_array[-1].append(sorted[i + 1])

          p +=1
     s += 1

outputFile = open ('summary.dist', 'w')
outputFile.write('min: ' + repr(global_min) + '\n')
outputFile.write('max: ' + repr(global_max) + '\n')
outputFile.write('cell: ' + repr(lattice_constant) + '\n')
outputFile.close()

# histogramin #

min_value = global_min
max_value = global_max
number_of_neighbours = int(number_of_particles - 1)

histogram = numpy.zeros((int(number_of_neighbours), int(number_of_bins + 1)), dtype=numpy.int)

bin_size = (max_value - min_value)/float(number_of_bins)

outputFile = open ('summary.hist', 'w')
outputFile.write('min used: ' + repr(min_value) + '\n')
outputFile.write('max used: ' + repr(max_value) + '\n')
outputFile.write('bin_size: ' + repr(bin_size) +  '\n')
outputFile.close()

counter = 0

length_distance_array = len(distance_array)

for neighbour_list in distance_array:
     for neighbour_index in range(number_of_neighbours):
          bin = int((neighbour_list[neighbour_index] - min_value)/bin_size)
          histogram[neighbour_index][bin] += 1
     counter += 1
#     print str(100*float(counter)/length_distance_array) + '%'

outputFile_nn_average = open ('nn_average.hist', 'w')
outputFile_nn_average.write('# bin (Angst), nn1, nn2, ... \n')


for bin_index in range(number_of_bins):

     outputFile_nn_average.write(repr( (bin_index*bin_size + bin_size/2.0 + min_value)*0.529177 ) + '  ')
     for neighbour_index in range(number_of_neighbours):
          outputFile_nn_average.write(str(histogram[neighbour_index][bin_index]/float(length_distance_array)) + ' ')
     outputFile_nn_average.write('\n')

outputFile_nn_average.close()
