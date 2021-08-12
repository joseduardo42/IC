from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve

#initial params
w = 100
h = 2
Ra = 10**3
C1 = 10*10**(-12)
C2 = 10**(-6)
RL = 50
Ra = 10**3
A = 100

F = np.array([[1, sin(2*pi*0*1/(2*h+1)),cos((2*pi*0*1/(2*h+1))), sin(2*pi*0*2/(2*h+1)), cos((2*pi*0*2/(2*h+1)))],
            [1, sin(2*pi*1*1/(2*h+1)),cos((2*pi*1*1/(2*h+1))), sin(2*pi*1*2/(2*h+1)), cos((2*pi*1*2/(2*h+1)))],
            [1, sin(2*pi*2*1/(2*h+1)),cos((2*pi*2*1/(2*h+1))), sin(2*pi*2*2/(2*h+1)), cos((2*pi*2*2/(2*h+1)))],
            [1, sin(2*pi*3*1/(2*h+1)),cos((2*pi*3*1/(2*h+1))), sin(2*pi*3*2/(2*h+1)), cos((2*pi*3*2/(2*h+1)))],
            [1, sin(2*pi*4*1/(2*h+1)),cos((2*pi*4*1/(2*h+1))), sin(2*pi*4*2/(2*h+1)), cos((2*pi*4*2/(2*h+1)))]])

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
        F_inv@((0.1*np.sign(non_linear))/((1 + (1.8/abs(non_linear))**5)**(1/5))) - (1/Ra)*Vb +C1*(omega@(Vb-Vc)),
        C2*omega@(Vb - Vc) -( 1/RL)*Vc 
    ])

#starting estimate and solve the system of nonlinear equations
amplitudes_guess = np.zeros(15)
y = fsolve(circuit_equations, amplitudes_guess)

print (y)
