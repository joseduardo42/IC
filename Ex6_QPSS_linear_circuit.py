import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from numpy import arange, concatenate, dot, linalg, pi
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
T1 = (1/f1)
T = (1/f2)

for h in range(2,10):

  #frequency -> time
  gamma_inv = np.array([[1] + [f(2*pi*(i)*(j+1)/(2*h+1)) for j in range(h) for f in (sin, cos)] for i in range(2*h+1)])

  #frequency -> time, one period ahead
  gamma_inv_T1 = np.array([[1] + [f(2*pi*f2*(j+1)*((i)*T/(2*h+1) + T1)) for j in range(h) for f in (sin, cos)] for i in range(2*h+1)])

  #time -> frequency
  gamma = linalg.inv(gamma_inv)

  #delay matrix
  D = gamma_inv_T1@gamma

  #vectors to storage results
  results_vc = []

  transient_result = np.zeros(2*h+1)

  def QPSS(V_shooting):

    shooting_voltage = np.zeros(2*h+1)
    for j in range(2*h+1):
        shooting_voltage[j] = V_shooting[j]
    #print (shooting_voltage)

  ################## shooting #####################
  ################## transient (5) #####################

    for i in range(2*h+1):

      t_sim = np.arange(i*T1, (i+1)*T1+deltat, deltat)

      if (t_sim[-1] != T1):
        t_sim = np.arange(i*T1, (i+1)*T1, deltat)

      #analysis in t0 = 0, T1, ..., (2N+1)T1

      vc0 = shooting_voltage[i]

      vs = Vm1*np.sin(2*pi*f1*t_sim[0])

      A1 = np.array([[R1, -R1],
                    [-R1, (R1 + R2 + R3)]], np.float64)

      b = np.array([[vs], [-vc0]], np.float64)
      
      lin_SM_t0 = linalg.solve(A1,b)

      i0 = float (lin_SM_t0[1]) #current i(0), i(T1),..., (2N)T1

      for t in np.delete(t_sim, 0):
        
          vs = Vm1*np.sin(2*pi*f1*t) #voltage source
          
          #linear system to solve in each t  
          A2 = np.array([[R1, -R1, 0],
              [-R1, (R1 + R2 + R3), 1],
              [0, -deltat/(2*C), 1]], np.float64)
          
          b = np.array([[vs], [0], [vc0 + i0*deltat/(2*C)]], np.float64)

          z = linalg.solve(A2, b) #solution of system in z

          #capacitor current 
          i0 = float (z[1])

          #capacitor voltage
          vc0 = float (z[2])

      transient_result[i] = vc0

    return np.concatenate([
      transient_result - D@shooting_voltage
      ])

  #solve QPSS function
  amplitudes_guess = np.zeros(2*h+1)
  y = fsolve(QPSS, amplitudes_guess)
  #print (y)

  t_sim = np.arange(0, 5*1/f2 + 1/(100 * f2), 1/(100 * f2))

  for t in t_sim:

      sinandcos = np.array([1] + [f(w2*(j+1)*t) for j in range(h) for f in (sin, cos)])
      Vc_time = dot(y, sinandcos)
      #print(Vc_time)

      results_vc.append (Vc_time)
  

  plt.plot (t_sim, results_vc, label = f'H = {h}')
  plt.title ('Tensão no capacitor')
  plt.ylabel ('(V)')
  plt.xlabel ('Tempo (milisegundos)')
  plt.grid()
  plt.legend()
  

#plot results
plt.show()