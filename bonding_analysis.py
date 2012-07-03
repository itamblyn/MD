#! /usr/bin/env python

import sys
import numpy

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print 'usage: ' + sys.argv[0] + ' TRAJEC.xyz [acell (Bohr)]'

def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
          if (input_value > 0): i+= 1
          if (input_value < 0): i-= 1
     return i

xyzFile = open(sys.argv[1], 'r')
number_of_particles = int(xyzFile.readline())

CELL_array = []
POSITION_array = []

for line in xyzFile:
    if len(line.split()) == 4:

         position_bonding_vector = [float(line.split()[1])/(0.529177), float(line.split()[2])/(0.529177),float(line.split()[3])/(0.529177), -10]
         POSITION_array.append(position_bonding_vector)

    elif len(line.split()) == 3:
         CELL_array.append([])
         for i in line.split():
             CELL_array[-1].append(float(i)/0.529177)
xyzFile.close()


number_of_snapshots = int(len(POSITION_array)/float(number_of_particles))

if len(sys.argv) == 3:
    acell = float(sys.argv[2])
    CELL_array = acell*numpy.ones((number_of_snapshots, 3), dtype=numpy.float)

sanity = (CELL_array[0][0]**2 + CELL_array[0][1]**2 + CELL_array[0][2]**2)**(1.0/2.0)

POSITION_array = numpy.array(POSITION_array)

dissFile = open('diss.dat', 'w')

s = 0 # counts over snapshots

while s < number_of_snapshots:

     molecule_number = 0
   
     p = 0 # counts over particles
     
     while p < number_of_particles:
          
          if (POSITION_array[s*number_of_particles + p][3] == -10):
     
               o = 0 # counts over other particles

               minimum_distance = sanity
          
               closest_to_p = p

               while o < number_of_particles:

                    dx = POSITION_array[s*number_of_particles + p][0] - POSITION_array[s*number_of_particles + o][0]
                    dy = POSITION_array[s*number_of_particles + p][1] - POSITION_array[s*number_of_particles + o][1]
                    dz = POSITION_array[s*number_of_particles + p][2] - POSITION_array[s*number_of_particles + o][2]

                    dx -= CELL_array[s][0]*pbc_round(dx/CELL_array[s][0])
                    dy -= CELL_array[s][1]*pbc_round(dy/CELL_array[s][1])
                    dz -= CELL_array[s][2]*pbc_round(dz/CELL_array[s][2])

                    distance = (dx**2 + dy**2 + dz**2)**(0.5)

                    if (distance > sanity): print 'Warning, problem with pbc'
               
                    if (distance < minimum_distance and distance != 0.0):
                         minimum_distance = distance
                         closest_to_p = o 
           
                    o +=1

               minimum_distance = sanity

               o = closest_to_p
               oo = 0

               if (POSITION_array[s*number_of_particles + o][3] == -10):
          
                    while oo < number_of_particles:
               
                         dx = POSITION_array[s*number_of_particles + o][0] - POSITION_array[s*number_of_particles + oo][0]
                         dy = POSITION_array[s*number_of_particles + o][1] - POSITION_array[s*number_of_particles + oo][1]
                         dz = POSITION_array[s*number_of_particles + o][2] - POSITION_array[s*number_of_particles + oo][2]

                         dx -= CELL_array[s][0]*pbc_round(dx/CELL_array[s][0])
                         dy -= CELL_array[s][1]*pbc_round(dy/CELL_array[s][1])
                         dz -= CELL_array[s][2]*pbc_round(dz/CELL_array[s][2])

                         distance = (dx**2 + dy**2 + dz**2)**(0.5)
                         
                         if (distance > sanity): print 'Warning, problem with pbc'
               
                         if (distance < minimum_distance and distance != 0.0):
                              minimum_distance = distance
                              closest_to_o = oo 
          
                         oo +=1

                    if (closest_to_p == o and closest_to_o == p):
                         molecule_number += 1
                         POSITION_array[s*number_of_particles + p][3] = o
                         POSITION_array[s*number_of_particles + o][3] = p
                         
                    else: POSITION_array[s*number_of_particles + p][3] = -1
                    
               else: POSITION_array[s*number_of_particles + p][3] = -1                    

          p +=1

     dissFile.write(str((2.0*float(molecule_number)/float(number_of_particles))) + '\n')

     s += 1

dissFile.close()

 ####################################
##                                  ##
## bonding analysis is now complete ##
##                                  ##
 ####################################

outputFile = open('TRAJEC.cbn', 'w')

outputFile.write('# comment = I like candy \n')
outputFile.write('# a = ' + str(acell) + '\n')
outputFile.write('# b = ' + str(acell) + '\n')
outputFile.write('# c = ' + str(acell) + '\n')
outputFile.write('# number_of_particles = ' + str(number_of_particles) + '\n')
outputFile.write('# number_of_neighbours = ' + str(1) + '\n')
outputFile.write('#\n')
outputFile.write('#\n')
outputFile.write('#\n')
outputFile.write('# units = bohr\n')

for s in range(number_of_snapshots):
    for p in range(number_of_particles):
        outputFile.write('Li ')
        outputFile.write(str(POSITION_array[s*number_of_particles + p][0]) + '  ')
        outputFile.write(str(POSITION_array[s*number_of_particles + p][1]) + '  ')
        outputFile.write(str(POSITION_array[s*number_of_particles + p][2]) + '  ')
        outputFile.write(str(int(POSITION_array[s*number_of_particles + p][3])) + '\n')
outputFile.close()
