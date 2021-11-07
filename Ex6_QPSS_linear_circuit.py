import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from numpy import linalg, pi
from scipy.optimize.minpack import fsolve

#ICs
R1 = 1
R2 = 2
R3 = 3
C = 1
Vm1 = 60
Vm2 = 40
f1 = 100
w1 = 2*pi*f1
f2 = 8
w2 = 2*pi*f2
deltat = 1/(100 * f1)
tf = (1/f1)
h = 2

#vectors to storage results
results_vc = []
shooting_voltage = np.zeros(5)

omega = np.zeros((5,5))
omega[1, 2] = -w2; omega[2, 1] = w2; omega[3, 4] = -2*w2; omega[4, 3] = 2*w2

#Define the system equations
def QPSS(V_shooting):

################## frequency #####################  
  def frequency_method(V):

    #vector of unknowns
    Va = np.zeros(5)
    Vb = np.zeros(5)
    Vc = np.zeros(5)
    for i in range(5):
      Va[i] = V[i]
      Vb[i] = V[5+i]
      Vc[i] = V[10+i]

    #definition of amplitude source and Va in time-domain
    A_amplitude = np.array([0, Vm2, 0, 0, 0])
    Va_aux = Va - A_amplitude
    vc1 = Vb - Vc

    return np.concatenate([
        Va_aux/R1 + (Va_aux - Vb)/R2,
        (Vb - Va_aux)/R2 - C*omega@vc1,
        C*omega@vc1 - Vc/R3
    ])
      
  amplitudes_guess = np.zeros(15)
  frequency_voltage = np.array(fsolve(frequency_method, amplitudes_guess))

################## shooting #####################
  for i in range(2*h+1):

    t_sim = np.arange(i*(tf+deltat), (i+1)*tf+deltat, deltat)
    #analysis in t0 = 0, T1, ..., (2N+1)T1

    #vc0 unknown of shooting function (5 unknowns)
    vc0 = V_shooting[0+i]
    
    #adjustment of vc0 for each interaction
    if (i > 0):
      vc0 = vc0f
      
    vs = Vm1*np.sin(2*pi*f1*float(t_sim[0]))

    A1 = np.array([[R1, -R1],
                    [-R1, (R1 + R2 + R3)]], np.float64)

    b = np.array([[vs], [-vc0]], np.float64)
    
    lin_SM_t0 = linalg.solve(A1,b)

    shooting_voltage[0+i] = vc0

    i0 = float (lin_SM_t0[1]) #current i(0)
    i = 1
    
    for t in np.delete(t_sim, 0):
        vs = Vm1*np.sin(2*pi*f1*t) #voltage source
        
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

        vc0 = vc
        vc0f = vc

        i+=1    

    print (vc0f)
  return np.concatenate([
    shooting_voltage - (frequency_voltage[5:10] - frequency_voltage[10:15])
    ])

#solve QPSS function
amplitudes_guess = np.zeros(5)
y = fsolve(QPSS, amplitudes_guess)

#waveforms of QPSS
t_sim = np.arange(0, 1/f2+1/(100 * f2), deltat)
for t in 2*t_sim:
    Vc_time =  y[0] +  y[1]*sin(w2*t)  + y[2]*cos(w2*t)  + y[3]*sin(2*w2*t)  + y[4]*cos(2*w2*t) 
    
    results_vc.append (Vc_time)

#plot results
plt.plot (t_sim, results_vc)
plt.title ('Tensão no capacitor')
plt.ylabel ('(V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()