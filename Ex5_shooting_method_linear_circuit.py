import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, linalg, pi
from scipy.optimize.minpack import fsolve

#ICs
R1 = 1
R2 = 2
R3 = 3
C = 10**-2
A = 220
f = 60
deltat = 1/(20 * f)
tf = 1/f
t0 = 0
t_sim = np.arange(0, tf+deltat, deltat) #vector to transient in shooting method

def shooting_method(x):
    
    #Analysis in t0
    vs = (A*np.sin(2*pi*f*t0))
    A1 = np.array([[R1, -R1],
                    [-R1, (R1 + R2 + R3)]], np.float64)
    b = np.array([[vs], [-x[0]]], np.float64)
    lin_MT_t0 = linalg.solve(A1,b)

    vc0 = float(x[0])
    i0 = float (lin_MT_t0[1]) #current i(0)
    
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
        vc0 = vc
        
    #discretized equation
    return [vc - x[0] - 1/(2*C)*deltat*i0- 1/(2*C)*deltat*ic]

y = fsolve(shooting_method, 0)

#initial condition to voltage in capacitor from the shooting method
vc0 = float(y[0])

t_plot = np.arange(0, tf+deltat, deltat) #vector with positions to plot
t_sim = np.arange(deltat, tf+deltat, deltat) #vector with time to simulate, without the t0 = 0

#vectors to store results
result_vc = []
result_ic = []
result_is = []

##Analysis in t0 (to plot)
A1 = np.array([[R1, -R1],
            [-R1, (R1 + R2 + R3)]], np.float64)
b = np.array([[0], [-vc0]], np.float64)

x = linalg.solve(A1,b)

i0 = float (x[1]) #current in t0 (to plot)

result_vc.append (vc0)
result_ic.append (i0)

#transient analysis equals activity 2
for t in t_sim:

    vs = (A*np.sin(2*pi*f*t))
 
    A2 = np.array([[R1, -R1, 0],
        [-R1, (R1 + R2 + R3), 1],
        [0, -deltat/(2*C), 1]], np.float64)
    b = np.array([[vs], [0], [vc0 + i0*deltat/(2*C)]], np.float64)
    
    x = linalg.solve(A2, b)

    i_s = float (x[0])
    result_is.append (i_s)

    ic = float (x[1])
    result_ic.append (ic)
    i0 = ic

    vc = float (x[2])
    result_vc.append (vc)
    vc0 = vc
        
plt.plot (t_plot, result_vc)
plt.title ('Tensão no capacitor')
plt.ylabel ('Tensão no capacitor (V)')
plt.xlabel ('Tempo (segundos)')
plt.grid()
plt.show ()
