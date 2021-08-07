from math import sin, cos, pi
import numpy as np
from numpy.core.fromnumeric import shape
from numpy.linalg import inv, solve
from scipy.optimize import fsolve

w = 100
h = 2
Ra = 10**3
C1 = 10*10**(-12)
C2 = 10**(-6)

F = np.array([[1, sin(2*pi*0*1/(2*h+1)),cos((2*pi*0*1/(2*h+1))), sin(2*pi*0*2/(2*h+1)), cos((2*pi*0*2/(2*h+1)))],
             [1, sin(2*pi*1*1/(2*h+1)),cos((2*pi*1*1/(2*h+1))), sin(2*pi*1*2/(2*h+1)), cos((2*pi*1*2/(2*h+1)))],
             [1, sin(2*pi*2*1/(2*h+1)),cos((2*pi*2*1/(2*h+1))), sin(2*pi*2*2/(2*h+1)), cos((2*pi*2*2/(2*h+1)))],
             [1, sin(2*pi*3*1/(2*h+1)),cos((2*pi*3*1/(2*h+1))), sin(2*pi*3*2/(2*h+1)), cos((2*pi*3*2/(2*h+1)))],
             [1, sin(2*pi*4*1/(2*h+1)),cos((2*pi*4*1/(2*h+1))), sin(2*pi*4*2/(2*h+1)), cos((2*pi*4*2/(2*h+1)))]])


def calc(IF):
    vetor_amplitudes_IF = np.array([IF[0], IF[1], IF[2],  IF[3], IF[4]])
    IFaux = F@vetor_amplitudes_IF
    
    return[((0.1*np.sign(IFaux)) / ((1 + (1.8/abs(IFaux))**5)**(1/5)))]

y = fsolve(calc, [0, 0, 0, 0, 0])
    """
def calc(V):
    return[
        V[0],
        V[1],
        V[2] - A,
        V[3],
        V[4],
        -V[10]/RL,
        C2*(w*V[12] - w*V[7]) - V[11]/RL,
        C2*(w*V[6] - w*V[11]) - V[12]/RL,
        C2*(2*w*V[14] - 2*w*V[9]) - V[13]/RL,
        C2*(2*w*V[8] - 2*w*V[13]) - V[14]/RL
    ]

y = fsolve(calc, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
"""
