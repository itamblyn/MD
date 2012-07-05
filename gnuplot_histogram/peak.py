#! /usr/bin/env python

import numpy 
import sys

print 'usage: ' + sys.argv[0] + ' nn_average.hist'

inputFile = open (sys.argv[1],'r')
inputFile.readline() # ignores the comment line at the top of the file...

histogram_array = []

#read line into array 
for line in inputFile.readlines():

    # add a new sublist
    histogram_array.append([])

    # loop over the elemets, split by whitespace
    for i in line.split():
    
        # convert to integer and append to the last
        # element of the list
        histogram_array[-1].append(float(i))

inputFile.close()

outputFile_peak_expectation_value = open ('peak_expectation_value.dat', 'w')
outputFile_peak_expectation_value_transpose = open ('peak_expectation_value_transpose.dat', 'w')
outputFile_peak_max_value = open ('peak_max_value.dat', 'w')
outputFile_peak_max_value_transpose = open ('peak_max_value_transpose.dat', 'w')
column = 1

while column < len(histogram_array[0]):
     row = 0
     expectation_value = 0
     max_value = 0
     while row < len(histogram_array):
          expectation_value += histogram_array[row][0]*histogram_array[row][column]
          if histogram_array[row][column] > max_value: 
               max_value = histogram_array[row][column]
               max_index = histogram_array[row][0]
          row +=1
     outputFile_peak_expectation_value.write(str(expectation_value) + '\n')     
     outputFile_peak_expectation_value_transpose.write(str(expectation_value) + ' ')
     outputFile_peak_max_value.write(str(max_index) + '\n')
     outputFile_peak_max_value_transpose.write(str(max_index) + ' ')
     
     column +=1

outputFile_peak_expectation_value_transpose.write('\n') 
outputFile_peak_max_value_transpose.write('\n')

outputFile_peak_expectation_value.close()
outputFile_peak_expectation_value_transpose.close()
outputFile_peak_max_value.close()
outputFile_peak_max_value_transpose.close()
