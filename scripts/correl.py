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
   
file1=sys.argv[1]
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

import numpy
print numpy.corrcoef(list1, list2)[0, 1]

