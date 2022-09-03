from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from Ex8_PAC_QPSS_nonlinear import h1, k, x1_1tom, w1

"""
This code have as objective the harmonic balance analysis of a non linear
circuit. To calculate the harmonic balance, it is necessary to provide the initial
conditions in the components of circuits and also the equations of circuit.
The simulation of transient depends of the wave amplitudes obtained in the
harmonic balance
"""

# params
C1 = 10 * 10 ** -12
C2 = 1 * 10 ** -6
R1 = 1 * 10 ** 3
RL = 50
f2 = 1.01 * 10 ** 9
w2 = 2 * pi * f2
Vm2 = 1
h2 = 1
k2 = 2 * (2 * h2 + 1)
T = (1 / f2)

# frequency -> time (F = gamma_inv)
gamma_inv = np.array([[1] + [f(2 * pi * i * (j + 1) / k) for j in range(h1) for f in (sin, cos)] for i in
                      range(k)])

# time -> frequency (F⁻¹ = gamma)
gamma = inv(gamma_inv)

# omega_2tons matrix
omega_2tons = np.zeros((k2, k2))
position_omega_vector = 0
for h in np.arange(-h2, h2 + 1, 1):
    omega_2tons[2 * position_omega_vector + 1, 2 * position_omega_vector] = - (w2 + h * w1)
    omega_2tons[2 * position_omega_vector, 2 * position_omega_vector + 1] = (w2 + h * w1)
    position_omega_vector += 1

# creating G_1tom (first matrix in right hand of (22))
non_linear = gamma_inv @ x1_1tom

g_1tom = gamma @ ((0.1 * np.sign(non_linear)) / ((1 + (1.8 / abs(non_linear)) ** 5) ** (1 / 5)))
for p in np.delete(range(len(g_1tom)), 0):
    g_1tom[p] = g_1tom[p]/2

G_1tom = np.zeros((k2, k2))
for i in range(k2):
    if i % 2 == 0:
        r1 = range(i, -1, -1)
        r2 = range(i + 2, k2)
    else:
        r1 = range(i, k2)
        r2 = range(i-2, -1, -1)

    for (n, j) in enumerate(r1):
        G_1tom[i, j] = g_1tom[n]

    for (n, j) in enumerate(r2):
        if n % 2 == 0:
            G_1tom[i, j] = g_1tom[n+2]
        else:
            G_1tom[i, j] = -g_1tom[n]


def hb_lin(v):
    # vector of unknowns
    Va = np.zeros(k2)
    Vb = np.zeros(k2)
    Vc = np.zeros(k2)
    for j in range(k2):
        Va[j] = v[j]
        Vb[j] = v[k2 + j]
        Vc[j] = v[2 * k2 + j]

    # definition of amplitude source and Va in time-domain
    A_amplitude = np.zeros(k2)
    A_amplitude[2] = Vm2
    vc2 = Vb - Vc

    return np.concatenate([
        Va - A_amplitude,
        G_1tom @ Va - (1 / R1) * Vb
        - (C1 * omega_2tons) @ Vb - (C2 * omega_2tons) @ vc2,
        C2 * omega_2tons @ vc2 - Vc / RL
    ])


amplitudes_guess = np.zeros(3 * k2)
linear_result = fsolve(hb_lin, amplitudes_guess)
x1_lin = linear_result[:k2]
x2_lin = linear_result[k2: 2 * k2]
x3_lin = linear_result[2 * k2: 3 * k2]
X_c1_lin = (C1 * omega_2tons) @ x2_lin
