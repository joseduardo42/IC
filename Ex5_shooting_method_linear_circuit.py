import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg, pi
from scipy.optimize.minpack import fsolve

#ICs
R1 = 1
R2 = 2
R3 = 3
C = 1
A = 100
f = 15.9155
deltat = 1/(20 * f)
tf = (1/f)
t0 = 0

#time domain
t_sim = np.arange(deltat, tf + deltat, deltat) #vector to transient in shooting method
t_plot = np.arange(0, tf + deltat, deltat) #vector with positions to plot

#global vector
result_vc = np.zeros(len(t_plot))

#shooting method
def shooting_method(x):
    #Analysis in t0
    vs = A*np.sin(2*pi*f*t0)
    vc0 = x[0]

    A1 = np.array([[R1, -R1],
                    [-R1, (R1 + R2 + R3)]], np.float64)
    b = np.array([[vs], [-vc0]], np.float64)

    lin_SM_t0 = linalg.solve(A1,b)

    result_vc[0] = vc0
    i0 = float (lin_SM_t0[1]) #current i(0)

    i = 1
    #transient
    for t in t_sim:
        vs = A*np.sin(2*pi*f*t) #voltage source
        
        #linear system to solve in each t  
        A2 = np.array([[R1, -R1, 0],
            [-R1, (R1 + R2 + R3), 1],
            [0, -deltat/(2*C), 1]], np.float64)
            
        b = np.array([[vs], [0], [vc0 + i0*deltat/(2*C)]], np.float64)

        z = linalg.solve(A2, b) #solution of system in z

        #capacitor current 
        ic = float (z[1])
        i0 = ic
        #capacitor voltage
        vc = float (z[2])
        result_vc[i] = vc
        vc0 = vc
        i+=1

    #discretized equation
    return [result_vc[-1] - x[0]]

#shooting method
sm_result = fsolve(shooting_method, 30)
print (sm_result)

#plotting 1 cycle of transient
plt.plot (t_plot, result_vc)
plt.title ('Tensão no capacitor')
plt.ylabel ('Tensão no capacitor (V)')
plt.xlabel ('Tempo (segundos)')
plt.grid()
plt.show ()
