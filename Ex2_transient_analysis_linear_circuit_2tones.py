import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg, pi
 
"""
This code have as objective the transient analysis of a linear circuit with
two tones source.To calculate the transient analysis, it is necessary to 
provide the initialconditions in the components of circuits and also the 
equations of circuit.The simulation depends on previous conditions in each
new interation.
"""

#ICs
R1 = 1
R2 = 2
R3 = 3
C = 1
Vm1 = 60
Vm2 = 40
f1 = 100
f2 = 10
deltat = 1/(100 * f1)
tf = 3*(1/f2)
vc0 = 0.17842857

t_sim = np.arange(0, tf+deltat, deltat) #time vector to simulation, without t0
if (t_sim[-1] != tf):
    t_sim = np.arange(0, tf, deltat)

vs = Vm1*np.sin(2*pi*f1*t_sim[0]) + Vm2*np.sin(2*pi*f2*t_sim[0])#voltage in source in actual time (t=0)

#The matrix from MNA
A1 = np.array([[1, 0, 0, 0, 0],
               [(1/R1)+(1/R2), -1/R2, 0, 1, 0],
               [-1/R2, 1/R2, 0, 0, -1],
               [0, 0, -1/R3, 0, 1],
               [0, 1, -1, 0, 0]])
          
b = np.array([[vs], [0], [0], [0], [vc0]], np.float64)

x = linalg.solve(A1, b)#solving the system of linear equations in t x = [Va, Vb, Vc, is, ic1]

i0 = float (x[4]) #current in t0, to use in interations in
print (i0)

result_vc = [] #vector to storage the voltage at capacitor
result_ic = [] #vector to storage the current at capacitor
result_is = [] #vector to storage the current at source

#storage the parameters of cicuits at t0
result_vc.append (vc0)
result_ic.append (i0)
result_is.append (float (x[3]))
print (float(x[3]))

for t in np.delete(t_sim, 0):

    vs = Vm1*np.sin(2*pi*f1*t) + Vm2*np.sin(2*pi*f2*t)#voltage in source in actual time

    #The matrix from MNA         
    A2 = np.array([[1, 0, 0, 0, 0],
                   [(1/R1)+(1/R2), -1/R2, 0, 1, 0],
                   [-1/R2, 1/R2, 0, 0, -1],
                   [0, 0, -1/R3, 0, 1],
                   [0, 1, -1, 0, -deltat/(2*C)]], np.float64)

    b = np.array([[vs], [0], [0], [0], [vc0 + i0*deltat/(2*C)]], np.float64)

    x = linalg.solve(A2, b)  #solving the system of linear equations in t x = [Va, Vb, Vc, is, ic1]
    
    #current at source
    i_s = float (x[3])
    result_is.append (i_s)
    
    #current at capacitor
    ic = float (x[4])
    result_ic.append (ic)
    i0 = ic
    
    #voltage at capacitor
    vc = float (x[1] - x[2])
    result_vc.append (vc)
    vc0 = vc
    #print (i_s)

#plot transient analysis
plt.plot (t_sim, result_vc)
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
