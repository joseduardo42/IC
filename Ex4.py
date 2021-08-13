from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

f = 15.91549
deltat = 1/(20 * f)
tf = 5*1/f

#initial params
w = 100
h = 2
Ra = 10**3
C1 = 10*10**(-12)
C2 = 10**(-6)
RL = 50
Ra = 10**3
A = 100

F = np.array([[1, sin(2*pi*0*1/(2*h+1)), cos((2*pi*0*1/(2*h+1))), sin(2*pi*0*2/(2*h+1)), cos((2*pi*0*2/(2*h+1)))],
              [1, sin(2*pi*1*1/(2*h+1)), cos((2*pi*1*1/(2*h+1))), sin(2*pi*1*2/(2*h+1)), cos((2*pi*1*2/(2*h+1)))],
              [1, sin(2*pi*2*1/(2*h+1)), cos((2*pi*2*1/(2*h+1))), sin(2*pi*2*2/(2*h+1)), cos((2*pi*2*2/(2*h+1)))],
              [1, sin(2*pi*3*1/(2*h+1)), cos((2*pi*3*1/(2*h+1))), sin(2*pi*3*2/(2*h+1)), cos((2*pi*3*2/(2*h+1)))],
              [1, sin(2*pi*4*1/(2*h+1)), cos((2*pi*4*1/(2*h+1))), sin(2*pi*4*2/(2*h+1)), cos((2*pi*4*2/(2*h+1)))]])

F_inv = inv(F)

#omega matrix
omega = np.zeros((5,5))
omega[1, 2] = -w; omega[2, 1] = w; omega[3, 4] = -2*w; omega[4, 3] = 2*w

#Define the system equations
def circuit_equations(V):

    #vector of unknowns
    Va = np.zeros(5)
    Vb = np.zeros(5)
    Vc = np.zeros(5)
    for i in range(5):
        Va[i] = V[i]
        Vb[i] = V[5+i]
        Vc[i] = V[10+i]

    #definition of amplitude source and Va in time-domain
    A_amplitude = np.array([0, 0, A, 0, 0])
    non_linear = F@Va

    return np.concatenate([
        Va - A_amplitude,
        F_inv@((0.1*np.sign(non_linear))/((1 + (1.8/abs(non_linear))**5)**(1/5))) - (1/Ra)*Vb + C1*(omega@(Vb-Vc)),
        C2*omega@(Vb - Vc) - (1/RL)*Vc 
    ])

#starting estimate and solve the system of nonlinear equations
amplitudes_guess = np.zeros(15)
y = fsolve(circuit_equations, amplitudes_guess)
print (y)

#transient analysis
t_sim = np.arange(0, tf+deltat, deltat)#tempo de simulação
#vetores para armazenar valores de tensão
results_va = []
results_vb = []
results_vc2 = []

#cálculo das formas de onda do HB
for t in t_sim:
    Va_tempo =  y[0] + y[1]*sin(w*t) + y[2]*np.cos(w*t) + y[3]*sin(2*w*t) + y[4]*cos(2*w*t) 
    Vb_tempo = y[5] + y[6]*sin(w*t) + y[7]*cos(w*t) + y[8]*sin(2*w*t) + y[9]*cos(2*w*t) 
    Vc_tempo = y[10] + y[11]*sin(w*t) + y[12]*cos(w*t) + y[13]*sin(2*w*t) + y[14]*cos(2*w*t) 

    Vc2_tempo = (y[5] - y[10]) + (y[6]- y[11])*sin(w*t) + (y[7] - y[12])*cos(w*t) + (y[8] - y[13])*sin(2*w*t) + (y[9] - y[14])*cos(2*w*t) 
    
    results_va.append (Va_tempo)
    results_vb.append (Vb_tempo)
    results_vc2.append (Vc2_tempo)

plt.plot (t_sim, results_va)
plt.title ('Tensão nó A')
plt.ylabel ('(V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()
plt.plot (t_sim, results_vb)
plt.title ('Tensão no Capacitor 1')
plt.ylabel ('(V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()

plt.plot (t_sim, results_vc2)
plt.title ('Tensão no Capacitor 2')
plt.ylabel ('(V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()
