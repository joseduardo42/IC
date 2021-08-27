import numpy as np
from numpy import linalg
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from math import pi

#ICs
Ra = 10**3
RL = 50
C1 = 10**(-11)
C2 = 10**(-6)
#Ativ1 source
f = 15.9155
A = 100
deltat = 1/(20 * f)
tf = (1/f)

#time domain
t_sim = np.arange(deltat, tf + deltat, deltat)
t_plot = np.arange(0, tf + deltat, deltat) 
result_vc1 = np.zeros(len(t_plot))
result_vc2 = np.zeros(len(t_plot))
result_ic1 = np.zeros(len(t_plot))
result_ic2 = np.zeros(len(t_plot))

def shooting_method_non_linear(z):
    print (result_vc1,result_vc2 )
    #Analysis in t0
    Vc10 = float (z[0])
    Vc20 = float (z[1])

    A1 = np.array([[Ra, 0],
                    [0, RL]], np.float64)
    b = np.array([[-Vc10], [Vc10 - Vc20]])
    x = linalg.solve(A1, b)
  
    ic10 = float (x[0] - x[1])#current C1(0)
    ic20 = float (x[1])#current C2 i(0)

    result_vc1[0] = Vc10
    result_vc2[0] = Vc20
    result_ic1[0] = ic10
    result_ic2[0] = ic20
    #variables for interactions
    y1 = float (x[0])
    y2 = float (x[1])
    #transient
    i = 1
    for t in t_sim:

        Vs = A*np.cos(2*pi*f*t)

        #non-linear system to solve in each t
        def func(y):
            return[Ra*(y[0]-((0.1*np.sign(Vs)) / ((1 + (1.8/abs(Vs))**5)**(1/5)))) + Vc10 + (1/(2*C1))*deltat*((y[0]-y[1]) + ic10),
                   -Vc10 - 1/(2*C1)*deltat*((y[0]-y[1]) + ic10) + Vc20 + 1/(2*C2)*deltat*(y[1] + ic20) + RL*y[1]]

        y = fsolve(func, [y1, y2])

        #valtage in capacitor for the next iteraction
        Vc10 = Vc10 + 1/(2*C1)*deltat*((y[0] - y[1]) + ic10)
        Vc20 = Vc20 + 1/(2*C2)*deltat*(y[1] + ic20)

        #initial guess for the next interaction
        y1 = float (y[0])
        y2 = float (y[1])
        
        #variables for next interaction
        ic10 = y1 - y2
        ic20 = y2

        #capacitors voltage 
        result_vc1[i] = Vc10
        result_vc2[i] = Vc20
        i+=1

    #discretized equations
    return[result_vc1[-1] - z[0],
           result_vc2[-1] - z[1]]

#solving shooting method
sm_nonlin_result = fsolve(shooting_method_non_linear,np.array([100, 100]))

#plotting 1 cycle of transient in C1
plt.plot (t_plot, result_vc2)
plt.title ('Tensão no capacitor 1')
plt.ylabel ('Tensão no capacitor (V)')
plt.xlabel ('Tempo (segundos)')
plt.grid()
plt.show ()

#plotting 1 cycle of transient in C2
#plt.plot (t_plot, result_vc2)
#plt.title ('Tensão no capacitor')
#plt.ylabel ('Tensão no capacitor (V)')
#plt.xlabel ('Tempo (segundos)')
#plt.grid()
#plt.show ()
