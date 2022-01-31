from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

from Ex3_transien_nonlinear_circuit import nonlinear_element as nonlinear_element_transient
# from Ex7_PAC_HB_superposition import Vm1, f1, T1, f2, w1, k, h, R1, C1, C2, RL
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
Vm2 = 3
h = 2
k = 2 * h + 1
"""
This code have as objective the harmonic balance analysis of a non linear
circuit. To calculate the harmonic balance, it is necessary to provide the initial
conditions in the components of circuits and also the equations of circuit.
The simulation of transient depends of the wave amplitudes obtained in the
harmonic balance
"""

# frequency -> time (F = gamma_inv)
gamma_inv = np.array([[1] + [f(2 * pi * i * (j + 1) / k) for j in range(h) for f in (sin, cos)] for i in
                      range(k)])
# time -> frequency (F⁻¹ = gamma)
gamma = inv(gamma_inv)

# omega matrix
omega = np.zeros((k, k))
for i in range(h):
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
amplitudes_guess = np.zeros(3*k)
y = fsolve(circuit_equations, amplitudes_guess)

X_va = y[:k]
X_vb = y[k: 2 * k]
X_vc = y[2 * k: 3 * k]

n = int(100 * f1 / f2)
(t_sim, deltat) = np.linspace(0, 10 * (1 / f1), n, retstep=True)

results_va = []
results_vb = []
results_vc = []
nonlinear_element = []

# waveforms of HB
for t in t_sim:
    sinandcos = np.array([1] + [f(w1 * (j + 1) * t) for j in range(h) for f in (sin, cos)])
    Va_time = sinandcos @ X_va
    Vb_time = sinandcos @ X_vb
    Vc_time = sinandcos @ X_vc
    dependent_source = (0.1 * np.sign(Va_time)) / ((1 + (1.8 / abs(Va_time)) ** 5) ** (1 / 5))

    results_va.append(Va_time)
    results_vb.append(Vb_time)
    results_vc.append(Vc_time)
    nonlinear_element.append(dependent_source)

MSE = mean_squared_error(nonlinear_element, nonlinear_element_transient)
print(MSE)

# plot results
plt.plot(t_sim, nonlinear_element, label='HB')
plt.legend(loc="upper right")
plt.plot(t_sim, nonlinear_element_transient, label='Transitório')
plt.legend(loc="upper right")
plt.title('Corrente da fonte controlada')
plt.ylabel('(V)')
plt.xlabel('Tempo (mili segundos)')
plt.grid()
plt.show()
