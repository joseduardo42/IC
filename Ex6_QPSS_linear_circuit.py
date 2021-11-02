import numpy as np
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
deltat = 1/(20 * f1)
tf = (1/f1)
vs = 0
h = 2

F = np.array([[1, sin(2*pi*0*1/(2*h+1)), cos((2*pi*0*1/(2*h+1))), sin(2*pi*0*2/(2*h+1)), cos((2*pi*0*2)/(2*h+1))],
              [1, sin(2*pi*1*1/(2*h+1)), cos((2*pi*1*1/(2*h+1))), sin(2*pi*1*2/(2*h+1)), cos((2*pi*1*2)/(2*h+1))],
              [1, sin(2*pi*2*1/(2*h+1)), cos((2*pi*2*1/(2*h+1))), sin(2*pi*2*2/(2*h+1)), cos((2*pi*2*2)/(2*h+1))],
              [1, sin(2*pi*3*1/(2*h+1)), cos((2*pi*3*1/(2*h+1))), sin(2*pi*3*2/(2*h+1)), cos((2*pi*3*2)/(2*h+1))],
              [1, sin(2*pi*4*1/(2*h+1)), cos((2*pi*4*1/(2*h+1))), sin(2*pi*4*2/(2*h+1)), cos((2*pi*4*2)/(2*h+1))]])

t_sim = np.arange(0, tf+deltat, deltat)#time simulation
#vectors to storage results
results_va = []
results_vb = []
results_vc = []
capacitor_voltage = []
shooting_voltage = np.zeros(5)

omega = np.zeros((5,5))
omega[1, 2] = -w2; omega[2, 1] = w2; omega[3, 4] = -2*w2; omega[4, 3] = 2*w2

#Define the system equations
def QPSS(V_shooting):
  
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

    #analysis in t0

    #vc0 unknown of shooting function (5 unknowns)
    vc0 = V_shooting[0+i]
    vs = 0
    A1 = np.array([[R1, -R1],
                    [-R1, (R1 + R2 + R3)]], np.float64)

    b = np.array([[vs], [-vc0]], np.float64)
    
    lin_SM_t0 = linalg.solve(A1,b)

    shooting_voltage[0+i] = vc0

    capacitor_voltage.append (vc0)

    i0 = float (lin_SM_t0[1]) #current i(0)
    i = 1

    for t in t_sim:
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
        capacitor_voltage.append (vc)
        
        vc0 = vc
        i+=1

    #vc0 to 2N+1 transient
    vc0 = capacitor_voltage[-1]
  
  return np.concatenate([
    shooting_voltage - (frequency_voltage[5:10] - frequency_voltage[10:15])
    ])

amplitudes_guess = np.zeros(5)
x = fsolve(QPSS, amplitudes_guess)