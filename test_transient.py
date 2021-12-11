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

C = 1
I1 = 60
I2 = 40
f1 = 100
f2 = 10
deltat = 1/(100 * f1)
tf = (1/f2)
vc0 = -57.83404023

#The matrix from the mesh analysis in circuit

i_saux = I1*np.sin(2*pi*f1*0) + I2*np.sin(2*pi*f2*0) #current in t0, to use in interations in
print (i_saux)
t_sim = np.arange(deltat, tf+deltat, deltat) #time vector to simulation, without t0
t_plot_trans = np.arange(0, tf+deltat, deltat)  #vector to plot in each time of simulation

result_vc_trans = [] #vector to storage the voltage at capacitor
result_ic = [] #vector to storage the current at capacitor
result_is = [] #vector to storage the current at source

#storage the parameters of cicuits at t0
result_vc_trans.append (vc0)
result_ic.append (i_saux)
result_is.append (i_saux)

for t in t_sim:

        i_s = I1*np.sin(2*pi*f1*t) + I2*np.sin(2*pi*f2*t)#voltage in source in actual time      
        #system of mesh analysis to solve in actual time         
        vc = vc0 + (deltat/C)*((i_s + i_saux)/2)
        
        result_vc_trans.append (vc)
        vc0 = vc
        
        i_saux = i_s
        #print (i_s)

#plot transient analysis
plt.plot (t_plot_trans, result_vc_trans)
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
