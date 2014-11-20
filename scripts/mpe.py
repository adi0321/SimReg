'''
Computes Pearson Correlation
'''

from decimal import *
import sys
import os
import logging
import random
import math
import itertools
from operator import itemgetter, attrgetter


# scipy/numpy.
import numpy as np


### script ###    

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

if len(sys.argv) < 3:
    print 'Usage:\n\
    argv[1] - file 1 \n\
    argv[2] - file 2'
    sys.exit()
   
#true frequencies
file1=sys.argv[1]
#estimated frequencies
file2=sys.argv[2]

#import pandas
#Need to install pandas for this
	#colnames = ['names', 'freq']
	#data = pandas.read_csv(file1, names=colnames)


# Load File 1
list1 = np.loadtxt(file1) #load list1 from file into array
#print 'Print File 1:'
#print list1

list2 = np.loadtxt(file2) #load list2 from file into array

#Normalize all values
list1Norm = list1/np.sum(list1)
list2Norm = list2/np.sum(list2)

#print list1Norm
#print list2Norm
with np.errstate(divide='ignore', invalid='ignore'): 
    # division errors suppressed only within this block
    mpe = abs(list1Norm-list2Norm)/list1Norm
    mpe[list1Norm == 0] = 1000000
    #mpe[list1Norm == 0.00000000e+00] = 1000000

#print mpe

#mpeNorm = mpe/np.sum(mpe)

#print np.median(abs(list1Norm-list2Norm)/list1Norm)
#print np.median(mpeNorm)
print np.median(mpe)

