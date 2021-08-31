from typing import final
import numpy as np
from numpy import linalg
from scipy.optimize import fsolve, root
import matplotlib.pyplot as plt
from math import pi

#from Ex3 import t_plot__aux, result_vc1_transi,aux, result_vc2_transi

#ICs
Ra = 10**3
RL = 50
C1 = 10**(-11)
C2 = 10**(-6)
#Ativ1 source
f = 15.915
A = 100
deltat = 1/(100 * f)
tf = (1/f)
final_resnorm = 0

#time domain
t_sim = np.arange(deltat, tf + deltat, deltat)
t_plot = np.arange(0, tf + deltat, deltat)
result_vc1 = np.zeros(len(t_plot))
result_vc2 = np.zeros(len(t_plot))
result_ic1 = np.zeros(len(t_plot))
result_ic2 = np.zeros(len(t_plot))

#shooting method
def shooting_method_non_linear(z):
    final_resnorm = 0
    Vc10 = z[0]
    Vc20 = z[1]
    def func(x):
        return[ Ra*x[0] + Vc10,
                 -Vc10 + Vc20 + x[1]*RL]

    a = fsolve(func, [0, 0])
    #variables for interactions
    y1 = float (a[0])
    y2 = float (a[1])   

    ic10 = y1 - y2 #current C1 i(0)
    ic20 = y2 #current C2 i(0)

    result_vc1[0] = Vc10
    result_vc2[0] = Vc20
    result_ic1[0] = ic10
    result_ic2[0] = ic20

    #transient
    i = 1
    for t in t_sim:
        
        Vs = A*np.cos(2*pi*f*t)
        #non-linear system to solve in each t
        def func(y):

            return[Ra*(y[0]-((0.1*np.sign(Vs)) / ((1 + (1.8/abs(Vs))**5)**(1/5)))) + (Vc10 + (1/(2*C1))*deltat*((y[0] - y[1]) + ic10)),
                  -(Vc10 + 1/(2*C1)*deltat*((y[0] - y[1]) + ic10)) + (Vc20 + 1/(2*C2)*deltat*(y[1] + ic20)) + RL*y[1]]

        y = fsolve(func, [y1, y2], full_output = True)
        
        #resnorm internal fsolve
        resnorm = sum(y[1]['fvec']**2)
        if resnorm > final_resnorm:
            final_resnorm = resnorm

        #initial guess for the next interaction
        y1 = float (y[0][0])
        y2 = float (y[0][1])

        #valtage in capacitor for the next iteraction
        Vc10 = (Vc10 + 1/(2*C1)*deltat*((y1 - y2) + ic10))
        Vc20 = (Vc20 + 1/(2*C2)*deltat*(y2 + ic20))
        
        #variables for next interaction
        ic10 = y1 - y2
        ic20 = y2

        #capacitors voltage 
        result_vc1[i] = Vc10
        result_vc2[i] = Vc20
        result_ic1[i] = ic10
        result_ic2[i] = ic20
        i+=1

    print(final_resnorm)
    #discretized equations
    return[result_vc1[-1] - z[0],
           result_vc2[-1] - z[1]]

#solving shooting method
sm_nonlin_result = fsolve(shooting_method_non_linear, [100, 100], full_output= True)

#resnorm external fsolve

resnorm = sum(sm_nonlin_result[1]['fvec']**2)
if resnorm > final_resnorm:
    final_resnorm = resnorm

print(final_resnorm)

##plotting 1 cycle of transient in C1
#plt.plot (t_plot, result_vc1)
##plt.plot (t_plot__aux, result_vc1_transi[-aux:])
#plt.title ('Tensão no capacitor 1')
#plt.ylabel ('Tensão no capacitor 1 (V)')
#plt.xlabel ('Tempo (segundos)')
#plt.grid()
#plt.show ()
#
##plotting 1 cycle of transient in C2
#plt.plot (t_plot, result_vc2)
#plt.plot (t_plot__aux, result_vc2_transi[-aux:])
#plt.title ('Tensão no capacitor 2')
#plt.ylabel ('Tensão no capacitor (V)')
#plt.xlabel ('Tempo (segundos)')
#plt.grid ()
#plt.show ()
