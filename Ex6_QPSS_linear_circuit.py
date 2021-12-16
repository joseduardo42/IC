import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from numpy import linalg, pi
from scipy.optimize.minpack import fsolve
# from Ex2_transient_analysis_linear_circuit_2tones import result_vc_trans, t_sim_trans

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
  shooting_Va = np.zeros(n_unknows)
  shooting_Vb = np.zeros(n_unknows)
  shooting_Vc = np.zeros(n_unknows)
  shooting_is = np.zeros(n_unknows)
  shooting_ic = np.zeros(n_unknows)

  transient_Va = np.zeros(n_unknows)
  transient_Vb = np.zeros(n_unknows)
  transient_Vc = np.zeros(n_unknows)
  transient_is = np.zeros(n_unknows)
  transient_ic = np.zeros(n_unknows)

  def QPSS(V_shooting):

    #vector of unknowns
    for i in range(n_unknows):
        shooting_Va[i] = V_shooting[i]
        shooting_Vb[i] = V_shooting[n_unknows+i]
        shooting_Vc[i] = V_shooting[2*n_unknows+i]
        shooting_is[i] = V_shooting[3*n_unknows+i]
        shooting_ic[i] = V_shooting[4*n_unknows+i]

  ################## shooting #####################
  ################## transient (5) #####################

    for i in range(n_unknows):
      
      t_sim = np.arange(i*(T/(n_unknows)), i*(T/(n_unknows)) + T1 + deltat, deltat)
      if (t_sim[-1] != i*(T/(n_unknows)) + T1):
        
        t_sim = np.arange(i*(T/(n_unknows)), i*(T/(n_unknows)) + T1, deltat)

      ################## first logic ###################

      vs = Vm1*np.sin(2*pi*f1*t_sim[0]) + Vm2*np.sin(2*pi*f2*t_sim[0])
      
      vc0 = float (shooting_Vb[i] - shooting_Vc[i])
      A1 = np.array([[1, 0, 0, 0, 0],
                      [(1/R1)+(1/R2), -1/R2, 0, 1, 0],
                      [-1/R2, 1/R2, 0, 0, 1],
                      [0, 0, -1/R3, 0, 1],
                      [0, 1, -1, 0, 0]], np.float64)

      b = np.array([[vs], [0], [0], [0], [vc0]], np.float64)      
      
      z = linalg.solve(A1, b) #solution of system in z
      i0 = float (z[4])
      print(i0)
      
      for t in np.delete(t_sim, 0):
        
        vs = Vm1*np.sin(2*pi*f1*t) + Vm2*np.sin(2*pi*f2*t)#voltage source
        
        #linear system to solve in each t  
        A2 = np.array([[1, 0, 0, 0, 0],
                      [(1/R1)+(1/R2), -1/R2, 0, 1, 0],
                      [-1/R2, 1/R2, 0, 0, 1],
                      [0, 0, -1/R3, 0, 1],
                      [0, 1, -1, 0, -deltat/(2*C)]], np.float64)

        b = np.array([[vs], [0], [0], [0], [vc0 + i0*deltat/(2*C)]], np.float64)

        z = linalg.solve(A2, b) #solution of system in z

        #node voltages
        Va_t = float (z[0])
        Vb_t = float (z[1])
        Vc_t = float (z[2])
        #source current
        i_s = float (z[3])   
        #capacitor current 
        i0 = float (z[4])
        #capacitor voltage
        vc0 = Vb_t - Vc_t

      #transient results
      transient_Va[i] = Va_t
      transient_Vb[i] = Vb_t
      transient_Vc[i] = Vc_t
      transient_is[i] = i_s
      transient_ic[i] = i0

    return np.concatenate([
      transient_Va - D@shooting_Va,
      transient_Vb - D@shooting_Vb,
      transient_Vc - D@shooting_Vc,
      transient_is - D@shooting_is,
      transient_ic - D@shooting_ic
      ])
  
  #solve QPSS function
  amplitudes_guess = np.zeros(n_unknows*5)
  y = fsolve(QPSS, amplitudes_guess, full_output=True)

  #resnorm external fsolve
  resnorm = sum(y[1]['fvec']**2)
  if resnorm > final_resnorm:
      final_resnorm = resnorm
  #print(final_resnorm)
  # print (final_resnorm)

Y = gamma@(y[0][n_unknows:2*n_unknows]-y[0][2*n_unknows:3*n_unknows])

t_sim = np.array( [i*T1 for i in range (int(f1/f2) + 1)] )
# print(Y)
results_vc = []
for t in t_sim:

    sinandcos = np.array([1] + [f(w2*(j+1)*t) for j in range(h) for f in (sin, cos)])
    Vc_time = Y@sinandcos
    # print(Vc_time)
    results_vc.append (Vc_time)

#plot results
# plt.plot (t_sim, results_vc, "ob")
# plt.plot (t_sim_trans, result_vc_trans)
# plt.title ('Tensão no capacitor')
# plt.ylabel ('(V)')
# plt.xlabel ('Tempo (milisegundos)')
# plt.grid()
# plt.show()