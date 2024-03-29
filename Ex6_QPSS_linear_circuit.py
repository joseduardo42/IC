import numpy as np
from numpy import linalg, pi
import matplotlib.pyplot as plt
from math import sin, cos, pi
from scipy.optimize.minpack import fsolve
from Ex2_transient_analysis_linear_circuit_2tones import t_sim_trans as t_plot, result_vc_trans

# ICs
R1 = 1
R2 = 2
R3 = 3
C = 1
Vm1 = 60
Vm2 = 40
f1 = 100
w1 = 2 * pi * f1
f2 = 8
w2 = 2 * pi * f2
final_resnorm = 0

T1 = (1 / f1)
T = (1 / f2)

for h in range(4, 5):

    n_unknows = 2 * h + 1
    # frequency -> time
    gamma_inv = np.array([[1] + [f(2 * pi * i * (j + 1) / (2 * h + 1)) for j in range(h) for f in (sin, cos)] for i in
                          range(2 * h + 1)])

    # frequency -> time, one period ahead
    gamma_inv_T1 = np.array(
        [[1] + [f(2 * pi * f2 * (j + 1) * (i * T / (2 * h + 1) + T1)) for j in range(h) for f in (sin, cos)] for i in
         range(2 * h + 1)])

    # time -> frequency
    gamma = linalg.inv(gamma_inv)

    # delay matrix
    D = gamma_inv_T1 @ gamma

    # vectors to storage results
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


    def qpss(shooting_tension):

        # vector of unknowns

        for i in range(n_unknows):
            shooting_Va[i] = shooting_tension[i]
            shooting_Vb[i] = shooting_tension[n_unknows + i]
            shooting_Vc[i] = shooting_tension[2 * n_unknows + i]
            shooting_is[i] = shooting_tension[3 * n_unknows + i]
            shooting_ic[i] = shooting_tension[4 * n_unknows + i]

        ################## shooting #####################
        ################## transient (5) #####################

        for i in range(n_unknows):
            n = int(10 * f1 / f2)
            (t_sim, deltat) = np.linspace(i * (T / (2 * h + 1)), i * (T / (2 * h + 1)) + T1, n, retstep=True)

            vs = Vm1 * np.sin(2 * pi * f1 * t_sim[0]) + Vm2 * np.sin(2 * pi * f2 * t_sim[0])

            vc0 = shooting_Vb[i] - shooting_Vc[i]

            A1 = np.array([[1, 0, 0, 0, 0],
                           [(1 / R1) + (1 / R2), -1 / R2, 0, 1, 0],
                           [-1 / R2, 1 / R2, 0, 0, 1],
                           [0, 0, -1 / R3, 0, 1],
                           [0, 1, -1, 0, 0]], np.float64)

            b = np.array([[vs], [0], [0], [0], [vc0]], np.float64)

            z = linalg.solve(A1, b)  # solution of system in z

            i0 = float(z[4])

            global Va, Vb, Vc

            for t in np.delete(t_sim, 0):
                vs = Vm1 * np.sin(2 * pi * f1 * t) + Vm2 * np.sin(2 * pi * f2 * t)  # voltage source

                # linear system to solve in each t
                a2 = np.array([[1, 0, 0, 0, 0],
                               [(1 / R1) + (1 / R2), -1 / R2, 0, 1, 0],
                               [-1 / R2, 1 / R2, 0, 0, 1],
                               [0, 0, -1 / R3, 0, 1],
                               [0, 1, -1, 0, -deltat / (2 * C)]], np.float64)

                b = np.array([[vs], [0], [0], [0], [vc0 + i0 * deltat / (2 * C)]], np.float64)

                z = linalg.solve(a2, b)  # solution of system in z

                # node voltages
                Va = float(z[0])
                Vb = float(z[1])
                Vc = float(z[2])
                # source current
                i_s = float(z[3])
                # capacitor current
                i0 = float(z[4])
                # capacitor voltage
                vc0 = Vb - Vc

            # transient results
            transient_Va[i] = Va
            transient_Vb[i] = Vb
            transient_Vc[i] = Vc
            transient_is[i] = i_s
            transient_ic[i] = i0

        return np.concatenate([
            transient_Va - D @ shooting_Va,
            transient_Vb - D @ shooting_Vb,
            transient_Vc - D @ shooting_Vc,
            transient_is - D @ shooting_is,
            transient_ic - D @ shooting_ic
        ])


    # solve qpss function
    amplitudes_guess = np.zeros(n_unknows * 5)
    y = fsolve(qpss, amplitudes_guess, full_output=True)

    # resnorm external fsolve
    resnorm = sum(y[1]['fvec'] ** 2)
    if resnorm > final_resnorm:
        final_resnorm = resnorm

    Y = gamma @ (y[0][n_unknows:2 * n_unknows] - y[0][2 * n_unknows:3 * n_unknows])

    t_sim_waveform = np.array([i * T1 for i in range(int(f1 / f2) + 1)])

    results_vc = []

    for t_waveform in t_sim_waveform:
        sinandcos = np.array([1] + [f(w2 * (j + 1) * t_waveform) for j in range(h) for f in (sin, cos)])
        Vc_time = Y @ sinandcos
        results_vc.append(Vc_time)

    # plot results
    plt.plot(t_sim_waveform, results_vc, "ob", label=f'H = {h}')
    plt.plot(t_plot, result_vc_trans)
    plt.title('Tensão no capacitor')
    plt.ylabel('(V)')
    plt.xlabel('Tempo (segundos)')
    plt.grid()
    plt.legend()
plt.show()
