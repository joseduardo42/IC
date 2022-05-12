from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve

# params
C1 = 10 * 10 ** -12
C2 = 1 * 10 ** -6
R1 = 1 * 10 ** 3
RL = 50
f1 = 1 * 10 ** 9
w1 = 2 * pi * f1
f2 = 1.1 * 10 ** 9
w2 = 2 * pi * f2
Vm1 = 5
Vm2 = 0.2
h1 = 16
k = 2 * h1 + 1

# frequency -> time (F = gamma_inv)
gamma_inv = np.array([[1] + [f(2 * pi * i * (j + 1) / k) for j in range(h1) for f in (sin, cos)] for i in
                      range(k)])
# time -> frequency (F⁻¹ = gamma)
gamma = inv(gamma_inv)

# omega matrix
omega = np.zeros((k, k))
for i in range(h1):
    omega[2 * i + 1, 2 * i + 2] = - (i + 1) * w1
    omega[2 * i + 2, 2 * i + 1] = (i + 1) * w1


# Define the system equations
def circuit_equations(v):
    # vector of unknowns
    Va = np.zeros(k)
    Vb = np.zeros(k)
    Vc = np.zeros(k)

    for j in range(k):
        Va[j] = v[j]
        Vb[j] = v[k + j]
        Vc[j] = v[2 * k + j]

    # definition of amplitude source and Va in time-domain
    A_amplitude = np.zeros(k)
    A_amplitude[1] = Vm1
    non_linear = gamma_inv @ Va
    vc2 = Vb - Vc

    return np.concatenate([
        Va - A_amplitude,
        (gamma @ ((0.1 * np.sign(non_linear)) / ((1 + (1.8 / abs(non_linear)) ** 5) ** (1 / 5)))) - (1 / R1) * Vb
        - (C1 * omega) @ Vb - (C2 * omega) @ vc2,
        C2 * omega @ vc2 - Vc / RL
    ])


# starting estimate and solve the system of nonlinear equations
amplitudes_guess = np.zeros(3 * k)
nonlinear_result = fsolve(circuit_equations, amplitudes_guess)

x1_1tom = nonlinear_result[:k]
x2_1tom = nonlinear_result[k: 2 * k]
x3_1tom = nonlinear_result[2 * k: 3 * k]
X_c1_1tom = (C1 * omega) @ x2_1tom

