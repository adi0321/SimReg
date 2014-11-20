'''
uses clustering and least-squares to solve deconvolution
plot sum of absolute error (residuals)
'''
#Library used: http://openopt.org/CVXOPT
### imports ###

# system
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

# cvxpy
import cvxpy

import logging
log = logging.getLogger('qp')
log.setLevel(logging.NOTSET)
#log.setLevel(logging.ERROR)
#log.setLevel(logging.DEBUG)
#DEBUG print all logs
formatter = logging.Formatter('[%(levelname)s] %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
log.addHandler(handler)

### definitions ###

### classes ###

### internal functions ###

def _cqp(D, o, k, m, N):
    ''' solves using cvxpy '''
    
    # cast to object.
    D = cvxpy.matrix(D)
    N = cvxpy.matrix(N)	
    o = cvxpy.matrix(o)
    
    # create variables.
    f = cvxpy.variable(k, 1, name='f')
    x=cvxpy.variable(m, 1, name='x')
    y=cvxpy.variable(m, 1, name='y')    
	
	# create constraints.
    geqs = cvxpy.greater_equals(f,0.0)
	
	#TO DO: Sum of all f = sum of observed reads classes (and not equal to 1)
    sum1 = cvxpy.equals(cvxpy.sum(f), 1.0)
    #sum1 = cvxpy.equals(cvxpy.sum(f), sum_obs_freq)
    
    #3	
    #dev = cvxpy.equals(D*f-o-x,0.0)	
    
	#4. matrix N (m x m) * x - y = 0
    sizeConstr = cvxpy.equals(N*x-y,0.0)
    #Check now to do N^2
	#sizeConstr = cvxpy.equals(N^2*x-y,0.0)
    #This might not work but try
    #sizeConstr = cvxpy.equals(x/N-y,0.0)
	
    #constrs = [geqs, sum1, dev, sizeConstr]
    constrs = [geqs, sum1]
		
    log.debug('\tin _cqp function: \n\t\tPrint matrices shapes:')
    log.debug('\t\t\t%s', D.shape)
    log.debug('\t\t\t%s', f.shape)
    log.debug('\t\t\t%s', o.shape)
    
    # create the program.
    #p = cvxpy.program(cvxpy.minimize(cvxpy.norm2(y)),constraints=constrs)
    p = cvxpy.program(cvxpy.minimize(cvxpy.norm2(D*f-o)),constraints=constrs)
    p.options['abstol'] = 1e-6 ## 'abstol' - Absolute accuracy	Default: 1e-7
    p.options['reltol'] = 1e-5 ## 'reltol' - Relative accuracy	Default: 1e-6
    p.options['feastol'] = 1e-5 ## 'feastol' - Tolerance for feasibility conditions	Default: 1e-6
    p.options['maxiters'] = 500 ## 'maxiters' - Maximum number of iterations	Default: 100
    
    
    # solve the program.
    p.solve(quiet=True)
	
    # return results.
    #print np.around(f.value, decimals=20)
    
    #print "Print using loop"
    getcontext().prec = 20
    #for i in f.value:
    #    temp_fi=str(i).strip('[]')	
    #    print temp_fi
    
    return f.value

### script ###    

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

if len(sys.argv) < 6:
    print 'Usage:\n\
    argv[1] - number of transcripts \n\
    argv[2] - number of read classes \n\
    argv[3] - d_values file \n\
    argv[4] - o_values file \n\
	argv[5] - tr_lengths file \n\
	argv[6] - file with the size of each read class'
    sys.exit()
   


# dimensions.
k = int(sys.argv[1])
m = int(sys.argv[2])
d_values=sys.argv[3]
o_values=sys.argv[4]
tr_lengths_file=sys.argv[5]
rcSize_file=sys.argv[6]

#print ('# of transcripts k=' + str(k) + '\n and # read classes m=' + str(m) + '\n')

# get transcripts length
tr_length = np.loadtxt(tr_lengths_file) #load matrix from file into array

log.debug('Print transcripts length: %s', tr_length)
 
# get size of each read class
N = np.loadtxt(rcSize_file) #load matrix from file into array
log.debug('Print Read Classes Size: \n%s', N)

log.info('Load matrix D from file into array')
# D matrix
D = np.loadtxt(d_values) #load matrix from file into array
log.debug('Matrix D: \n%s', D)
 
# observed frequencies
o = np.zeros((m,1), dtype=np.float)
temp_o = np.loadtxt(o_values)
log.debug('temp_o: \n%s',temp_o)

log.debug('Print size of o_values file: \n%s', len(np.atleast_1d(temp_o)))
#o_values may be just one row

log.debug('Length of D: \n%s', len(np.atleast_1d(D)))

#This assumption was wrong because we can have only one o value but 2 transcripts
if (len(np.atleast_1d(temp_o)) == 1) and (len(np.atleast_1d(D)) == 1):
	#If single transcript --> Skip the solver and give 100% to it
	fileName, fileExtension = os.path.splitext(d_values)
	outfile=fileName+".iso.estimates"
	f_file = open( outfile, 'w' )
	f_file.write("[ 1.]\n")
	f_file.close()
	sys.exit(0)
elif len(np.atleast_1d(temp_o)) == 1 :
	o[0,0] = temp_o
else:
	for i in range(len(np.atleast_1d(temp_o))):
		o[i,0] = temp_o[i]

log.debug('Print observed read frequencies: \n%s', o)

#print '\nEstimate transcript frequencies:'
# estimate transcript frequencies
f = _cqp(D, o, k, m, N)
 
f=np.around(f, decimals=8)
#print np.around(f, decimals=6)
log.debug("Print f': \n%s\n-",f)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Print original f values (f') to file (before the values are divided by transcript length)
fileName2, fileExtension2 = os.path.splitext(d_values)
outfile2=fileName2+"_iso.estimates.original"
f_file2 = open( outfile2, 'w' )
f_file2.write("%s\n" % f)

f_file2.close()
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


################################################
#print "\nCompute sum of all f'"
sum_f_prime=0
for i in range(len(f)):
    sum_f_prime+=f[i]

#print ('Sum = ' + str(sum_f_prime) +'\n')

#################################################

#print "\n Compute sum of all (f' DIVIDED by tr_length)"

sum_f_prime_tr=0
getcontext().prec = 20
for i in range(len(f)):
    #print np.around(f[i]/tr_length[i], decimals=20)
    temp_fi=str(f[i]).strip('[]')
    #print temp_fi
    #print tr_length[i]	
    #print len(f)
    sum_f_prime_tr+=Decimal(temp_fi)/Decimal(tr_length[i])
#print "Sum = ",sum_f_prime_tr

#sys.exit()	

fileName, fileExtension = os.path.splitext(d_values)
outfile=fileName+".iso.estimates"
f_file = open( outfile, 'w' )

#print "Print final frequency: (f' / length) / sum_of_all "

for i in range(len(f)):
    temp_fi=str(f[i]).strip('[]')
    #print (Decimal(temp_fi)/Decimal(tr_length[i]))/Decimal(sum_f_prime_tr)
    iso_estimates=(Decimal(temp_fi)/Decimal(tr_length[i]))/Decimal(sum_f_prime_tr)
    f_file.write("%s\n" % iso_estimates)

f_file.close()
