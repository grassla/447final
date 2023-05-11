#### PHY447 Final Project by Alexis Grassl ####

#### This program calculates the Clebsch-Gordan Coefficients for a selection of quantum numbers where:
#### j1 is the total angular momentum for state 1,
#### j2 is the total angular momentum for state 2,
#### m1 is the magnetic quantum number for state 1,
#### m2 is the magnetic quantum number for state 2,
#### J is the combined total angular momentum of the states,
#### and M is the combined magnetic quantum number for the states

#### When choosing values for j1 and j2, this program will give all of the CGCs for every possible value of m1 and m2,
#### given that -j1<=m1<=j1 and -j2<=m2<=j2

#### NOTE that, if running in Jupyter notebook, it will not support interactive version or show 'Pretty Table.'

import math
import sympy
from sympy.physics.quantum.cg import CG
from sympy import S
from sympy import factorial
from sympy import sqrt
import numpy as np
import pandas as pd


#j1 = 1/2
#j2 = 3/2
#j1 = float(input('Enter a value for j1: '))
#j2 = float(input('Enter a value for j2: '))
#J = int(j1+j2)

while True:
    try:
        j1 = float(input('Enter a value for j1: '))
        j2 = float(input('Enter a value for j2: '))
        J = int(j1+j2)
        CGC = clebsch_gordan(j1, j2, m1, m2, J, M)
        break
    except:
        print('Invalid input. Try again.')
        continue
        
        
#print('  m1   ', ' m2   ', ' CGC ') #### Crude table ####
#print('----------------------')

data = []
for m1 in np.arange(-j1, j1+1, 1):
    for m2 in np.arange(-j2, j2+1, 1):
        M = int(m1 + m2)
        if(-J <= M <= J):
            CGC = clebsch_gordan(j1, j2, m1, m2, J, M)
            CGC_strip = '{:.4f}'.format(CGC).rstrip('0').rstrip('.')
            #print(('{}   {}   {:.4f}'.format(m1, m2, CGC))) #### Crude table ####
            #### Pretty table: ####

            data.append([m1, m2, CGC_strip])

#### DataFrame from the list of dictionaries ####
df = pd.DataFrame(data, columns=['m1', 'm2', 'CGC'])
df = df.set_index(['m1', 'm2', 'CGC'], drop=True)

#### SOME PYTHON APPLICATIONS MAY ONLY SUPPORT 'df' RATHER THAN 'display(df)' #####
#df
display(df)



def clebsch_gordan(j1, j2, m1, m2, J, M):
    if(2*j1 != math.floor(2*j1) or 2*j2 != math.floor(2*j2) or 2*m1 != math.floor(2*m1) or 2*m2 != math.floor(2*m2)):
        raise ValueError('Inputs must be integers or half-integers. Try again.')
    elif(j1-j2 <= J <= j1+j2 and m1+m2 == M):
        Numerator = sqrt((2 * J + 1) * factorial(J + j1 - j2) * factorial(J - j1 +j2) * factorial(j1 + j2 - J) * factorial(j1 - m1) * factorial(j1 + m1) * factorial (j2 - m2) * factorial (j2 + m2) * factorial(J + M) * factorial(J - M))
        Denominator = sqrt(factorial(j1 + j2 + J + 1))
        
        #### Summation Calcs ####
        sum = 0
        for k in range(0, J + M + 2):
            try:
                n = (((-1) ** k) / (factorial(k) * factorial(j1 + j2 - J - k) * factorial(j1 - m1 - k) * factorial(j2 + m2 - k) * factorial(J - j2 + m1 + k) * factorial(J - j1 - m2 + k)))
            #### Add the current term to the sum in a increments ####
                sum += n
            except: pass
        CGC = (Numerator / Denominator) * sum 
    else:
        raise ValueError('You have an illegitimate value. Try again.')
        CGC = 0
    return CGC

CGC = clebsch_gordan(j1, j2, m1, m2, J, M)


