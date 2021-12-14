import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from numpy import concatenate, dot, linalg, pi
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
f2 = 10
w2 = 2*pi*f2
final_resnorm = 0

deltat = 1/(100 * f1)
T1 = (1/f1)
T = (1/f2)

for h in range(1,2):

  n_unknows = 2*h+1
  
  #frequency -> time
  gamma_inv = np.array([[1] + [f(2*pi*(i)*(j+1)/(2*h+1)) for j in range(h) for f in (sin, cos)] for i in range(2*h+1)])

  #frequency -> time, one period ahead
  gamma_inv_T1 = np.array([[1] + [f(2*pi*f2*(j+1)*((i)*T/(2*h+1) + T1)) for j in range(h) for f in (sin, cos)] for i in range(2*h+1)])

  #time -> frequency
  gamma = linalg.inv(gamma_inv)

  #delay matrix
  D = gamma_inv_T1@gamma

  #vectors to storage results
  shooting_voltage = np.zeros(n_unknows)
  transient_result = np.zeros(n_unknows)

  def QPSS(V_shooting):

    Va = np.zeros(n_unknows)
    Vb = np.zeros(n_unknows)
    Vc = np.zeros(n_unknows)
    for i in range(n_unknows):
        Va[i] = V_shooting[i]
        Vb[i] = V_shooting[n_unknows+i]
        Vc[i] = V_shooting[n_unknows+i]

    n = n_unknows*3
  ################## shooting #####################
  ################## transient (5) #####################

    for i in range(n_unknows):
      print(i)

      t_sim = np.arange(i*(T/(n_unknows)), i*(T/(n_unknows)) + T1 + deltat, deltat)
      if (t_sim[-1] != i*(T/(n_unknows)) + T1):
        
        t_sim = np.arange(i*(T/(n_unknows)), i*(T/(n_unknows)) + T1, deltat)

      A2 = np.array([[1, 0, 0, 0, 0],
                    [-1/R1-1/R2, 1/R2, 0, 1, 0],
                    [-1/R2, 1/R2, 0, 0, -1],
                    [0, 0, -1/R3, 0, 1],
                    [0, 1, -1, 0, 0]], np.float64)
          
      b = np.array([[vs], [0], [0], [0], [(vc0)]], np.float64)
      z = linalg.solve(A2, b) #solution of system in z
      
      #capacitor current 
      i0 = float (z[3])

      for t in np.delete(t_sim, 0):
        
        vs = Vm1*np.sin(2*pi*f1*t) + Vm2*np.sin(2*pi*f2*t)#voltage source
        
        #linear system to solve in each t  
        A2 = np.array([[1, 0, 0, 0, 0],
                      [(-1/R1-1/R2), 1/R2, 0, 1, 0],
                      [-1/R2, 1/R2, 0, 0, -1],
                      [0, 0, -1/R3, 0, 1],
                      [0, 1, -1, 0, -deltat/(2*C)]], np.float64)
        
        b = np.array([[vs], [0], [0], [0], [vc0 + i0*deltat/(2*C)]], np.float64)
        z = linalg.solve(A2, b) #solution of system in z
        #capacitor current 
        i0 = float (z[3])
        #capacitor voltage
        vc0 = float (z[2])

      transient_result[i] = vc0
      #print(shooting_voltage)
      
    return np.concatenate([
      transient_result - D@shooting_voltage
      ])
  
  #solve QPSS function
  amplitudes_guess = np.ones(n_unknows)
  y = fsolve(QPSS, amplitudes_guess, full_output=True)
  #print (y)

  #resnorm external fsolve
  resnorm = sum(y[1]['fvec']**2)
  if resnorm > final_resnorm:
      final_resnorm = resnorm
  #print(final_resnorm)
  print (y[0])

Y = gamma@y[0]

t_sim = np.array( [i*T1 for i in range (int(f1/f2) + 1)] )
print(t_sim)
results_vc = []
for t in t_sim:

    sinandcos = np.array([1] + [f(w2*(j+1)*t) for j in range(h) for f in (sin, cos)])
    Vc_time = dot(Y, sinandcos)
    #print(Vc_time)
    results_vc.append (Vc_time)

plt.plot (t_sim, results_vc)
plt.title ('Tensão no capacitor')
plt.ylabel ('(V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
#plot results
plt.show()