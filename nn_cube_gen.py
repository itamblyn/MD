#! /usr/bin/env python

import sys
import commands
import numpy

print 'this code is working, apart from the fact that the atoms and voxels are not properly lined up...'
print 'usage: ' + sys.argv[0] + ' number_of_bins lattice_constant(bohr) xyz_file'

def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
          if (input_value > 0): i+= 1
          if (input_value < 0): i-= 1
     return i

number_of_bins, lattice_constant, input_filename, snapshot_stride = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], 1 # should be a user input

number_of_x_bins, number_of_y_bins, number_of_z_bins = number_of_bins, number_of_bins, number_of_bins

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

###### histgram stuff

histogram = numpy.zeros( (number_of_x_bins, number_of_y_bins, number_of_z_bins), dtype=numpy.int)

x_step = (size_of_box_x)/float(number_of_x_bins)
y_step = (size_of_box_y)/float(number_of_y_bins)
z_step = (size_of_box_z)/float(number_of_z_bins)

SNAPSHOT_array = numpy.zeros((number_of_particles, 3), dtype=numpy.float)
AVERAGE_array  = numpy.zeros((number_of_particles, 3), dtype=numpy.float)

xyzFile.seek(0)			# return to beginning of file

s = 0
while s < number_of_snapshots:

    xyzFile.readline()          # skip number of particles
    xyzFile.readline()          # skip MD step or cell parameters

    for p in range(number_of_particles):

        line = xyzFile.readline()

        x = float(line.split()[1])/0.529177
        y = float(line.split()[2])/0.529177
        z = float(line.split()[3])/0.529177

        SNAPSHOT_array[p][0] = x - (int(x/size_of_box_x+number_of_snapshots+0.5)-number_of_snapshots)*size_of_box_x
        SNAPSHOT_array[p][1] = y - (int(y/size_of_box_y+number_of_snapshots+0.5)-number_of_snapshots)*size_of_box_y
        SNAPSHOT_array[p][2] = z - (int(z/size_of_box_z+number_of_snapshots+0.5)-number_of_snapshots)*size_of_box_z

    # At this point, SNAPSHOT_array is good to go!!!

    for p in range(number_of_particles):        # counts over particles 

        AVERAGE_array[p][0] += SNAPSHOT_array[p][0]
        AVERAGE_array[p][1] += SNAPSHOT_array[p][1]
        AVERAGE_array[p][2] += SNAPSHOT_array[p][2]

        bin_x = int((SNAPSHOT_array[p][0] - size_of_box_x/2.0)/x_step)
        bin_y = int((SNAPSHOT_array[p][1] - size_of_box_y/2.0)/y_step)
        bin_z = int((SNAPSHOT_array[p][2] - size_of_box_z/2.0)/z_step)

        histogram[bin_x][bin_y][bin_z] += 1

    s += snapshot_stride

outputFile_nn_cube = open ('nn.cube', 'w')

outputFile_nn_cube.write('COMMENT\n')
outputFile_nn_cube.write('COMMENT\n')
outputFile_nn_cube.write(str(number_of_particles) + ' ' + str(-size_of_box_x/2.0) + ' ' + str(-size_of_box_y/2.0) + ' ' + str(-size_of_box_z/2.0) + '\n')

outputFile_nn_cube.write(str(number_of_x_bins) + ' ' + str(size_of_box_x/float(number_of_x_bins)) + ' 0.000000 ' + ' 0.000000 ' + '\n')
outputFile_nn_cube.write(str(number_of_y_bins) + ' ' + ' 0.000000 ' + str(size_of_box_y/float(number_of_y_bins)) + ' 0.000000 ' + '\n')
outputFile_nn_cube.write(str(number_of_z_bins) + ' ' + ' 0.000000 ' + ' 0.000000 ' + str(float(size_of_box_z/float(number_of_z_bins))) + '\n')

for p in range(number_of_particles):
    outputFile_nn_cube.write('1 0.000000 ')
    outputFile_nn_cube.write(str(AVERAGE_array[p][0]/number_of_snapshots) + ' ' + str(AVERAGE_array[p][1]/number_of_snapshots) + ' ' + str(AVERAGE_array[p][2]/number_of_snapshots) + '\n')


for ix in range(number_of_x_bins): 
   for iy in range(number_of_y_bins):
      for iz in range(number_of_z_bins):
         voxel_value = float(histogram[ix][iy][iz])/(x_step*y_step*z_step)
         outputFile_nn_cube.write('% .6E' % voxel_value + ' ')
         if (iz % 6 == 5):
            outputFile_nn_cube.write('\n')
      outputFile_nn_cube.write('\n');

outputFile_nn_cube.close()
