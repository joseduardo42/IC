import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg, pi

"""
This code have as objective the transient analysis of a linear circuit.
To calculate the transient analysis, it is necessary to provide the initial
conditions in the components of circuits and also the equations of circuit.
The simulation depends on previous conditions in each new interation

"""

#ICs
R1 = 1
R2 = 2
R3 = 3
C = 1
A = 100
f = 10**9
deltat = 1/(20 * f)
tf = 1/f
vs = 0
vc0 = -0.1983515

#The matrix from the mesh analysis in circuit
A1 = np.array([[R1, -R1],
            [-R1, (R1 + R2 + R3)]], np.float64)
b = np.array([[vs], [-vc0]], np.float64)

x = linalg.solve(A1,b)

i0 = float (x[1]) #current in t0, to use in interations in 

t_sim = np.arange(deltat, tf+deltat, deltat) #time vector to simulation, without t0
t_plot = np.arange(0, tf+deltat, deltat) #vector to plot in each time of simulation
result_vc_trans = [] #vector to storage the voltage in capacitor
result_ic = [] #vector to storage the current at capacitor
result_is = [] #vector to storage the current at source
#storage the parameters of cicuits at t0
result_vc_trans.append (vc0)
result_ic.append (i0)
result_is.append (float (x[0]))

for t in t_sim:

    vs = (100*np.sin(2*pi*f*t)) #voltage in source in actual time
    #system of mesh analysis to solve in actual time        
    A2 = np.array([[R1, -R1, 0],
        [-R1, (R1 + R2 + R3), 1],
        [0, -deltat/(2*C), 1]], np.float64)
    b = np.array([[vs], [0], [vc0 + i0*deltat/(2*C)]], np.float64)
    
    x = linalg.solve(A2, b) #solving the system of linear equations
    
    #current at source
    i_s = float (x[0])
    result_is.append (i_s)
    
    #current at capacitor
    ic = float (x[1])
    result_ic.append (ic)
    i0 = ic
    
    #voltage at capacitor
    vc = float (x[2])
    result_vc_trans.append (vc)
    vc0 = vc

#variables for plot comparing to shooting method. Otherwise, use t_plot and result_vc       
t_plot__aux = np.arange(0, 1/f +deltat, deltat)
aux = len(t_plot__aux)

plt.plot (t_plot, result_vc_trans)
plt.title ('Tensão no capacitor')
plt.ylabel ('Tensão no capacitor (V)')
plt.xlabel ('Tempo (segundos)')
plt.grid()
plt.show ()

""""
plt.plot (t_plot, result_ic)
plt.title ('Corrente no capacitor')
plt.ylabel ('Corrente no capacitor (A)')
plt.xlabel ('Tempo (segundos)')
plt.grid()
plt.show ()

plt.plot (t_plot, result_is)
plt.title ('Tensão da fonte')
plt.ylabel ('Tensão (V)')
plt.xlabel ('Tempo (segundos)')
plt.grid()
plt.show ()
"""
