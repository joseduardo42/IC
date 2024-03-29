from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

"""
Thia code have as objective the harmonic balance analysis of a non linear
circuit. To calculate the harmonic balance, it is necessary to provide the initial
conditions in the components of circuits and also the equations of circuit.
The simulation of transient depends of the wave amplitudes obtained in the
harmonic balance
"""

# initial params
f = 10 ** 9
deltat = 1 / (20 * f)
tf = 5 * 1 / f
w = 2 * pi * f
h = 2
Ra = 10 ** 3
C1 = 10 * 10 ** (-12)
C2 = 10 ** (-6)
RL = 50
A = 5

F = np.array([[1, sin(2 * pi * 0 * 1 / (2 * h + 1)), cos((2 * pi * 0 * 1 / (2 * h + 1))),
               sin(2 * pi * 0 * 2 / (2 * h + 1)), cos((2 * pi * 0 * 2) / (2 * h + 1))],
              [1, sin(2 * pi * 1 * 1 / (2 * h + 1)), cos((2 * pi * 1 * 1 / (2 * h + 1))),
               sin(2 * pi * 1 * 2 / (2 * h + 1)), cos((2 * pi * 1 * 2) / (2 * h + 1))],
              [1, sin(2 * pi * 2 * 1 / (2 * h + 1)), cos((2 * pi * 2 * 1 / (2 * h + 1))),
               sin(2 * pi * 2 * 2 / (2 * h + 1)), cos((2 * pi * 2 * 2) / (2 * h + 1))],
              [1, sin(2 * pi * 3 * 1 / (2 * h + 1)), cos((2 * pi * 3 * 1 / (2 * h + 1))),
               sin(2 * pi * 3 * 2 / (2 * h + 1)), cos((2 * pi * 3 * 2) / (2 * h + 1))],
              [1, sin(2 * pi * 4 * 1 / (2 * h + 1)), cos((2 * pi * 4 * 1 / (2 * h + 1))),
               sin(2 * pi * 4 * 2 / (2 * h + 1)), cos((2 * pi * 4 * 2) / (2 * h + 1))]])

F_inv = inv(F)

# omega matrix
omega = np.zeros((5, 5))
omega[1, 2] = -w
omega[2, 1] = w
omega[3, 4] = -2 * w
omega[4, 3] = 2 * w


# Define the system equations
def circuit_equations(V):
    # vector of unknowns
    Va = np.zeros(5)
    Vb = np.zeros(5)
    Vc = np.zeros(5)
    for i in range(5):
        Va[i] = V[i]
        Vb[i] = V[5 + i]
        Vc[i] = V[10 + i]

    # definition of amplitude source and Va in time-domain
    A_amplitude = np.array([0, A, 0, 0, 0])
    non_linear = F @ Va
    vc2 = Vb - Vc

    return np.concatenate([
        Va - A_amplitude,
        (F_inv @ ((0.1 * np.sign(non_linear)) / ((1 + (1.8 / abs(non_linear)) ** 5) ** (1 / 5)))) - (1 / Ra) * Vb
        - (C1 * omega) @ Vb - (C2 * omega) @ vc2,
        C2 * omega @ vc2 - Vc / RL
    ])


# starting estimate and solve the system of nonlinear equations
amplitudes_guess = np.zeros(15)
y = fsolve(circuit_equations, amplitudes_guess)
print(y)

# transient analysis

t_sim = np.arange(0, tf + deltat, deltat)  # time simulation
# vectors to storage results
results_va = []
results_vb = []
results_vc = []

# waveforms of HB
for t in t_sim:
    Va_time = y[0] + y[1] * sin(w * t) + y[2] * cos(w * t) + y[3] * sin(2 * w * t) + y[4] * cos(2 * w * t)
    Vb_time = y[5] + y[6] * sin(w * t) + y[7] * cos(w * t) + y[8] * sin(2 * w * t) + y[9] * cos(2 * w * t)
    Vc_time = y[10] + y[11] * sin(w * t) + y[12] * cos(w * t) + y[13] * sin(2 * w * t) + y[14] * cos(2 * w * t)
    dependent_source = ((0.1 * np.sign(Va_time)) / ((1 + (1.8 / abs(Va_time)) ** 5) ** (1 / 5)))

    results_va.append(Va_time)
    results_vb.append(Vb_time)
    results_vc.append(Vc_time)

# plot results
plt.plot(t_sim, results_vb)
plt.title('Corrente da fonte controlada')
plt.ylabel('(V)')
plt.xlabel('Tempo (mili segundos)')
plt.grid()
plt.show()

# plt.plot (t_sim, results_vb)
# plt.title ('Tensão no Capacitor 1')
# plt.ylabel ('(V)')
# plt.xlabel ('Tempo (mili segundos)')
# plt.grid()
# plt.show()

# plt.plot (t_sim, results_vc2)
# plt.title ('Tensão no Capacitor 2')
# plt.ylabel ('(V)')
# plt.xlabel ('Tempo (mili segundos)')
# plt.grid()
# plt.show ()
