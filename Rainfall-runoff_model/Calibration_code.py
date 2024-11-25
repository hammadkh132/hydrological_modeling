# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 02:28:22 2024

@author: Hammad
"""
#Script for optimization code

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.optimize as opt
from permetrics.regression import RegressionMetric
from math import e

c = 0.172
Tr = 5.96
A = 247000000 # Area in m2
delta_t = 24*3600


filename = "C:/Users/Hammad/Downloads/Documents/Problem_set.csv" # Have to change the file path

data = pd.read_csv(filename, sep=",", encoding="utf-8", header=0)

data['I(m/sec)'] = data['I(mm/day)'] * (1/(1000*24*3600)) # conversion fromm mm/day to m/sec

 
# creating the list for guess of parameters
# guess = [0.172,5.96]
guess = [0.2, 10.0]

# Creating bounds (make sure this is a list of tuples)
# bnds = [(0.1, 0.30), (1, 10.0)]  # Each parameter has its own bounds
bnds = [(0.2, 0.3), (0.5, 50.0)]


I = data['I(m/sec)'].values

### flow function
def fun(params):
    c= params[0]
    Tr = params[1]
    Q_o = 0
    Qcal = [0]
    for i in range(len(I)-1):
        #Q_o = c*A*(data['I(m/sec)'][i]) + K
        Q = c*A*I[i] + (Q_o-c*I[i])*(np.exp(-delta_t/(Tr*24*3600)))
        Q_successive = c*A*I[i + 1] + (Q - c*A*I[i + 1])*(np.exp(-delta_t/(Tr*24*3600)))
        Q_o = Q
        Qcal.append(Q_successive)
    Qcal = np.array(Qcal)
    
    # Ensure the lengths match
    if len(Qcal) != len(I):
        raise ValueError("Qcal and I must have the same length.")
    
    evaluator = RegressionMetric(I, Qcal)
    n = evaluator.nash_sutcliffe_efficiency()
    return -n
        

# running optimization 
result = opt.minimize(fun, guess, method='Nelder-Mead', bounds=bnds)

print("Optimization Result:", result)

optimized_nse = -result.fun
print("Optimized NSE:", optimized_nse)


    