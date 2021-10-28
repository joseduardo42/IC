import matplotlib.pyplot as plt
import numpy as np
from math import sin, cos, pi
from numpy import linalg, pi
from scipy.optimize.minpack import fsolve
from lib.hb_functions import  omega
 
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
w1 = 2*pi*f1
f2 = 8
w2 = 2*pi*f2
deltat = 1/(20 * f1)
tf = (1/f1)
vs = 0


t_sim = np.arange(0, tf+deltat, deltat)#time simulation
#vectors to storage results
results_va = []
results_vb = []
results_vc = []
results_vc1 = []
result_vc_frequency = np.zeros(len(t_sim))

omega = np.zeros((5,5))
omega[1, 2] = -w2; omega[2, 1] = w2; omega[3, 4] = -2*w2; omega[4, 3] = 2*w2
#Define the system equations
def QPSS(y):
    
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

      return np.concatenate([
          Va_aux/R1 + (Va_aux - Vb)/R2,
          (Vb - Va_aux)/R2 - C*omega@(Vb - Vc),
          C*omega@(Vb - Vc) - Vc/R3
      ])

  amplitudes_guess = np.zeros(15)

  y = fsolve(frequency_method, amplitudes_guess)
  
  for t in t_sim:
    Va_time =  y[0] +  y[1]*sin(w1*t)  + y[2]*cos(w1*t)  + y[3]*sin(2*w1*t)  + y[4]*cos(2*w1*t) 
    Vb_time =  y[5] +  y[6]*sin(w1*t)  + y[7]*cos(w1*t)  + y[8]*sin(2*w1*t)  + y[9]*cos(2*w1*t)
    Vc_time =  y[10] + y[11]*sin(w1*t) + y[12]*cos(w1*t) + y[13]*sin(2*w1*t) + y[14]*cos(2*w1*t)
    
  results_va.append (Va_time)
  results_vb.append (Vb_time)
  results_vc.append (Vc_time)
  results_vc1.append (Vb_time - Vc_time)

  result_vc_frequency = []
  vs = Vm1*np.sin(2*pi*f1*0)

  vc0 = y[0]
  A1 = np.array([[R1, -R1],
                  [-R1, (R1 + R2 + R3)]], np.float64)
  b = np.array([[vs], [-vc0]], np.float64)
  lin_SM_t0 = linalg.solve(A1,b)
  result_vc_frequency[0] = vc0
  i0 = float (lin_SM_t0[1]) #current i(0)

  i = 1
  #transient
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
      result_vc_frequency[i] = vc
      vc0 = vc
      i+=1

  return[results_vc1[-1] - result_vc_frequency[-1]]

x = fsolve(QPSS, 30)
