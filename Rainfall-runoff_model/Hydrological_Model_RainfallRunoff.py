# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 23:11:50 2024

@author: Hammad Khalid
"""
# importing modules and libraries
import csv
import numpy as np
import pandas as pd
import math
from math import e, sqrt
import  matplotlib.pyplot as plt
from permetrics.regression import RegressionMetric

##### watershed area 
A = 247*(10**6) # m2
### Time period between two observations
delta_t = 24*3600

##### defining paramters
c = 0.172 #runoff coefficient
Tr = 5.96*24*3600 #Resident Time (days)

##### Coding the model

##### Imprting the data file

filename = "C:/Users/Hammad/Downloads/Documents/Problem_set.csv" #change the input file as per desired.
data = pd.read_csv(filename, sep=",", encoding="utf-8", header=0)

data['I(m/sec)'] = data['I(mm/day)'] * (1/(1000*24*3600))

##### Outflow Calculation
t=0
Q_o = 0
Qcal = [0]
for i in range(len(data['I(m/sec)'])-1):
    #Q_o = c*A*(data['I(m/sec)'][i]) + K
    Q = c*A*(data['I(m/sec)'][i]) + ((Q_o-c*A*data['I(m/sec)'][i])*(e**-(delta_t/Tr)))
    Q_successive = c*A*(data['I(m/sec)'][i+1]) + (Q-c*A*data['I(m/sec)'][i+1])*e**-(delta_t/Tr)
    Q_o = Q
    print (Q_successive)
    Qcal = np.append(Qcal, Q_successive)
data['Qcalculated(m3/sec)'] = Qcal

##### Plotting of Qmeasured and Qcalculated

fig, ax = plt.subplots(figsize = (7,4))
plt.plot(data["Date"], data["Qmeasured(m3/sec)"], color='r', label='Q_measured', linestyle='--', marker='o')
plt.plot(data["Date"], data["Qcalculated(m3/sec)"], color='b', label='Q_calculated', marker='o')
plt.xlabel("Date")
plt.ylabel("Q (m3/sec)")
plt.title("HYDROGRAPH")
ax.tick_params(axis='x', rotation=90)
plt.legend()
plt.show()

#### Calculating Efficiency of the model

#### Creating list
Qmeasured = [data['Qmeasured(m3/sec)']] 
Qcalculated = [data['Qcalculated(m3/sec)']] 


#### Root Mean sqaure error
evaluator = RegressionMetric(Qmeasured, Qcalculated)
print("Root Mean Square Error:\n") 
print(evaluator.root_mean_squared_error())
        
### Nash-Sutcliffe Efficiency  
print("\nNash Sutcliffe Efficiency:\n")
print(evaluator.nash_sutcliffe_efficiency())
     
#### Normalized Nash-Sutcliffe Efficiency   
print("\nNormalized Nash Sutcliffe Efficiency:\n")
print(evaluator.normalized_nash_sutcliffe_efficiency())

#### KGE - Kling-Gupta Efficiency
print("\nKling-Gupta Efficiency:\n")
print(evaluator.kling_gupta_efficiency())

#### Component of KGE
### Pearson Coefficient 
r = data['Qmeasured(m3/sec)'].corr(data['Qcalculated(m3/sec)'])
print("Pearson Coefficient (r) : ",r)
### relative variability
alpha =  ((data['Qcalculated(m3/sec)'].std()/data['Qcalculated(m3/sec)'].mean())/
          (data['Qmeasured(m3/sec)'].std()/data['Qmeasured(m3/sec)'].mean()))
print("Relative Variability : ",alpha)
### bias
B =  data['Qcalculated(m3/sec)'].mean()/data['Qmeasured(m3/sec)'].mean()
print("Bias : ",B)

### Manual KGE calculation
KGE_manual = 1-sqrt((1-r)**2+(B-1)**2+(alpha-1)**2)


# Plotting Q-Q plot

# Qmeasured = [data['Qmeasured(m3/sec)']] 
# Qcalculated = [data['Qcalculated(m3/sec)']] 

# percs = np.linspace(0,100,21)
# qn_a = np.percentile(Qmeasured, percs)
# qn_b = np.percentile(Qcalculated, percs)

# plt.plot(qn_a,qn_b, ls="", marker="o")

# x = np.linspace(np.min((qn_a.min(),qn_b.min())), np.max((qn_a.max(),qn_b.max())))
# plt.plot(x,x, color="k", ls="--")

# plt.show()






   





