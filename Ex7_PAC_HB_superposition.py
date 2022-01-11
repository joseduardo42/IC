from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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
f1 = 1 * 10 ** 9
w1 = 2 * pi * f1
f2 = 1.1 * 10 ** 9
w2 = 2 * pi * f2
Vm1 = 5
Vm2 = 3

T1 = (1 / f1)
T = (1 / f2)

for h in range(1, 2):
    k = 2 * h + 1
    # frequency -> time
    F = np.array([[1] + [f(2 * pi * i * (j + 1) / (2 * h + 1)) for j in range(h) for f in (sin, cos)] for i in
                  range(2 * h + 1)])
    # time -> frequency
    F_inv = inv(F)

    # omega matrix
    omega = np.zeros((5, 5))
    omega[1, 2] = -w
    omega[2, 1] = w
    omega[3, 4] = -2 * w
    omega[4, 3] = 2 * w

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
            (F_inv @ ((0.1 * np.sign(non_linear)) / ((1 + (1.8 / abs(non_linear)) ** 5) ** (1 / 5)))) - (
                        1 / R1) * Vb - (
                    C1 * omega) @ Vb - (C2 * omega) @ vc2,
            C2 * omega @ vc2 - (1 / RL) * Vc
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
