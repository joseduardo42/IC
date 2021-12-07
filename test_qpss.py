import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from numpy import dot, linalg, pi
from scipy.optimize.minpack import fsolve

#ICs
R1 = 1
R2 = 2
R3 = 3
C = 1
I1 = 60
I2 = 40
f1 = 100
w1 = 2*pi*f1
f2 = 8
w2 = 2*pi*f2
final_resnorm = 0

#[f1, f2] = [f2, f1]
#[w1, w2] = [w2, w1]
#print (f1,f2)

deltat = 1/(100 * f1)
T1 = (1/f1)
T = (1/f2)

for h in range(7,8):

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
  shooting_voltage = np.zeros(2*h+1)
  transient_result = np.zeros(2*h+1)

  def QPSS(V_shooting):

  ################## shooting #####################
  ################## transient (5) #####################

    for i in range(2*h+1):

      t_sim = np.arange(i*T1, (i+1)*T1+deltat, deltat)

      if (t_sim[-1] != T1):
        t_sim = np.arange(i*T1, (i+1)*T1, deltat)

      #analysis in t0 = 0, T1, ..., (2N+1)T1

      vc0 = V_shooting[i]
      shooting_voltage[i] = vc0

      i_saux = I1*np.sin(2*pi*f1*t_sim[0]) + I2*np.sin(2*pi*f2*t_sim[0])
      
      for t in np.delete(t_sim, 0):
        
        i_s = I1*np.sin(2*pi*f1*t) + I2*np.sin(2*pi*f2*t)#voltage in source in actual time

        #system of mesh analysis to solve in actual time         
        vc = vc0 + (deltat/C)*((i_s + i_saux)/2)
        
        vc0 = vc

        i_saux = i_s
        #print (i_s)
      
      transient_result[i] = vc0

    return np.concatenate([
      transient_result - D@shooting_voltage
      ])

  #solve QPSS function
  amplitudes_guess = np.ones(2*h+1)
  y = fsolve(QPSS, amplitudes_guess, full_output=True)
  
  #resnorm external fsolve
  resnorm = sum(y[1]['fvec']**2)
  if resnorm > final_resnorm:
      final_resnorm = resnorm
  print(final_resnorm)
  #print (y[0])

  t_sim = np.arange(0, 3*1/f2 + 1/(100 * f2), 1/(100 * f2))

  for t in t_sim:

      sinandcos = np.array([1] + [f(w2*(j+1)*t) for j in range(h) for f in (sin, cos)])
      Vc_time = dot(y[0], sinandcos)
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