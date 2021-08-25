import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, linalg, pi
from scipy.optimize.minpack import fsolve

#ICs
R1 = 1
R2 = 2
R3 = 3
C = 10**(-2)
A = 100
f = 60
deltat = 1/(20 * f)
tf = 1/f
t0 = 0

t_sim = np.arange(deltat, tf+deltat, deltat) #vector to transient in shooting method
t_plot = np.arange(0, tf+deltat, deltat) #vector with positions to plot
#global vector
result_vc = []

def shooting_method(x):
    #Analysis in t0
    vs = (A*np.sin(2*pi*f*t0))
    A1 = np.array([[R1, -R1],
                    [-R1, (R1 + R2 + R3)]], np.float64)
    b = np.array([[vs], [-x[0]]], np.float64)
    lin_SM_t0 = linalg.solve(A1,b)

    vc0 = float(x[0])
    i0 = float (lin_SM_t0[1]) #current i(0)
    
    for t in t_sim:
        vs = (A*np.sin(2*pi*f*t)) #voltage source
        
        #linear system to solve in each t  
        A2 = np.array([[R1, -R1, 0],
            [-R1, (R1 + R2 + R3), 1],
            [0, -deltat/(2*C), 1]], np.float64)
            
        b = np.array([[vs], [0], [-vc0 + i0*deltat/(2*C)]], np.float64)

        z = linalg.solve(A2, b) #solution of system in z

        #capacitor current 
        ic = float (z[1])
        i0 = ic
        #capacitor voltage
        vc = float (z[2])
        result_vc.append (vc)
        vc0 = vc

    #discretized equation
    return [result_vc[-1] - x[0]]

#shooting method
sm_result = fsolve(shooting_method, [0])
print (sm_result)

#numbers of elements of each interation
elements_last_interation = len(t_plot)

#plotting 1 cycle of transient
plt.plot (t_plot, result_vc[-elements_last_interation:])
plt.title ('Tensão no capacitor')
plt.ylabel ('Tensão no capacitor (V)')
plt.xlabel ('Tempo (segundos)')
plt.grid()
plt.show ()
