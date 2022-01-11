from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from Ex7_PAC_HB_superposition import Vm1, f1, T1, f2, w1, k, h, R1, C1, C2, RL

"""
This code have as objective the harmonic balance analysis of a non linear
circuit. To calculate the harmonic balance, it is necessary to provide the initial
conditions in the components of circuits and also the equations of circuit.
The simulation of transient depends of the wave amplitudes obtained in the
harmonic balance
"""

# frequency -> time
F = np.array([[1] + [f(2 * pi * i * (j + 1) / k) for j in range(h) for f in (sin, cos)] for i in
              range(k)])
# time -> frequency
F_inv = inv(F)

# omega matrix
omega = np.zeros((k, k))
for i in range(h):
    omega[2 * i + 1, 2 * i + 2] = - (i + 1) * w1
    omega[2 * i + 2, 2 * i + 1] = (i + 1) * w1

n = int(20 * f1 / f2)
(t_sim, deltat) = np.linspace(0, T1, n, retstep=True)


# Define the system equations
def circuit_equations(v):
    # vector of unknowns
    Va = np.zeros(5)
    Vb = np.zeros(5)
    Vc = np.zeros(5)
    for i in range(5):
        Va[i] = v[i]
        Vb[i] = v[5 + i]
        Vc[i] = v[10 + i]

    # definition of amplitude source and Va in time-domain
    A_amplitude = np.array([0, A, 0, 0, 0])
    non_linear = F @ Va
    vc2 = Vb - Vc

    return np.concatenate([
        Va - A_amplitude,
        (F_inv @ ((0.1 * np.sign(non_linear)) / ((1 + (1.8 / abs(non_linear)) ** 5) ** (1 / 5)))) - (1 / R1) * Vb - (
                C1 * omega) @ Vb - (C2 * omega) @ vc2,
        C2 * omega @ vc2 - (1 / RL) * Vc
    ])


# starting estimate and solve the system of nonlinear equations
amplitudes_guess = np.zeros(15)
y = fsolve(circuit_equations, amplitudes_guess)
print(y)