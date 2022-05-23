from turtle import color
import numpy as np
from numpy import linalg, pi
import matplotlib.pyplot as plt
from math import sin, cos, pi
from regex import R
from scipy.optimize.minpack import fsolve
from Ex7_superposition import t_sim as t_plot, results_vb

# ICs
C1 = 10 * 10 ** -12
C2 = 1 * 10 ** -6
R1 = 1 * 10 ** 3
RL = 50
f1 = 1 * 10 ** 9
f2 = 1.1 * 10 ** 9
fA = (f1 + f2)
fB = (abs(f1 - f2))
w1 = 2 * pi * f1 
w2 = 2 * pi * f2
Vm1 = 5
Vm2 = 0.2
h1 = 2
final_resnorm = 0

T1 = (1 / f1)
T = (1 / f2)
T1_nl = (1 / fA)
T_nl = (1 / fB)

for h in range(2, 3):
    k = 2 * h + 1
    # frequency -> time
    gamma_inv = np.array([[1] + [f(2 * pi * i * (j + 1) / k) for j in range(h) for f in (sin, cos)] for i in
                          range(k)])

    # frequency -> time, one period ahead
    gamma_inv_T1 = np.array(
        [[1] + [f(2 * pi * fB * (j + 1) * (i * T_nl / k + T1_nl)) for j in range(h) for f in (sin, cos)] for i in
         range(k)])

    # time -> frequency
    gamma = linalg.inv(gamma_inv)

    # delay matrix
    D = gamma_inv_T1 @ gamma

    # vectors to storage results
    shooting_Va = np.zeros(k)
    shooting_Vb = np.zeros(k)
    shooting_Vc = np.zeros(k)
    shooting_ic1 = np.zeros(k)
    shooting_ic2 = np.zeros(k)

    transient_Va = np.zeros(k)
    transient_Vb = np.zeros(k)
    transient_Vc = np.zeros(k)
    transient_ic1 = np.zeros(k)
    transient_ic2 = np.zeros(k)


    def qpss(shooting_voltage):

        # vector of unknowns

        for i in range(k):
            shooting_Va[i] = shooting_voltage[i]
            shooting_Vb[i] = shooting_voltage[k + i]
            shooting_Vc[i] = shooting_voltage[2 * k + i]
            shooting_ic1[i] = shooting_voltage[3 * k + i]
            shooting_ic2[i] = shooting_voltage[4 * k + i]

        ################## shooting #####################
        ################## transient (5) #####################

        for i in range(k):

            deltat = T1_nl / k
            t_sim = np.arange(i * (T_nl / k), i * (T_nl / k) + T1_nl, deltat)

            vs = Vm1 * np.sin(2 * pi * f1 * t_sim[0]) + Vm2 * np.sin(2 * pi * f2 * t_sim[0])

            Vc10 = shooting_Vb[i]
            Vc20 = shooting_Vb[i] - shooting_Vc[i]
            shooting_Va[i] = vs

            # system of nodal analysis to solve in actual time
            def func(x):

                if shooting_Va[i] == 0:
                    i_fnl = 0
                else:
                    i_fnl = ((0.1 * np.sign(shooting_Va[i])) / ((1 + (1.8 / abs(shooting_Va[i])) ** 5) ** (1 / 5)))

                return [i_fnl * R1 - x[0] * R1 + Vc10,
                        -Vc10 + Vc20 + x[1] * RL]

            initial_conditions = fsolve(func, np.array(
                [0, 0]))  # solving the system of nonlinear equations in t

            Va = shooting_Va[i]
            Vb = shooting_Vb[i]
            Vc = shooting_Vc[i]
            ic10 = float(initial_conditions[0] - float(initial_conditions[1]))
            ic20 = float(initial_conditions[1])

            for t in np.delete(t_sim, 0):
                vs = Vm1 * np.sin(2 * pi * f1 * t) + Vm2 * np.sin(2 * pi * f2 * t)

                # system of nodal analysis to solve in actual time
                def func(x):

                    if x[0] == 0:
                        i_fnl = 0
                    else:
                        i_fnl = ((0.1 * np.sign(x[0])) / ((1 + (1.8 / abs(x[0])) ** 5) ** (1 / 5)))

                    return [x[0] - vs,
                            i_fnl - x[1] / R1 - (
                            (2 * C1 / deltat) * (x[1] - Vc10) - ic10) - (
                            (2 * C2 / deltat) * ((x[1] - x[2]) - Vc20) - ic20),
                            -((2 * C2 / deltat) * ((x[1] - x[2]) - Vc20) - ic20) + x[2] / RL]

                z = fsolve(func, np.array(
                    [0, 0, 0]))  # solving the system of nonlinear equations in t

                # node voltages
                Va = float(z[0])
                Vb = float(z[1])
                Vc = float(z[2])

                # capacitor current
                ic10 = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5))) - Vb / R1 - Vc / RL
                ic20 = Vc / RL

                # capacitor voltage
                Vc10 = Vb
                Vc20 = Vb - Vc

            
            # transient results
            transient_Va[i] = Va
            transient_Vb[i] = Vb
            transient_Vc[i] = Vc
            transient_ic1[i] = ic10
            transient_ic2[i] = ic20

        print(Va)
        return np.concatenate([
            
            transient_Va - D @ shooting_Va,
            transient_Vb - D @ shooting_Vb,
            transient_Vc - D @ shooting_Vc,
            transient_ic1 - D @ shooting_ic1,
            transient_ic2 - D @ shooting_ic2
        ])
        
        
    
    # solve qpss function
    amplitudes_guess = [transient_Va, transient_Vb, transient_Vc, transient_ic1, transient_ic2]
    
    y = fsolve(qpss, amplitudes_guess, full_output=True)

    # resnorm external fsolve
    resnorm = sum(y[1]['fvec'] ** 2)
    if resnorm > final_resnorm:
        final_resnorm = resnorm


    Y = gamma @ (y[0][k:2 * k])


    t_sim_QPSS = np.array([i * (T_nl / k) for i in range(k)])
    results_vc = []

    for t_waveform in t_sim_QPSS:
        sinandcos = np.array([1] + [f(w1 * (j + 1) * t_waveform) for j in range(h) for f in (sin, cos)])
        Vc_time = Y @ sinandcos
        results_vc.append(Vc_time)

    # plot results
    plt.plot(t_sim_QPSS, transient_Va, "ob",color='red' , label=f'Chutes iniciais Va')
    plt.plot(t_sim_QPSS, transient_Vb, "ob",color='yellow' , label=f'Chutes iniciais Vb')
    plt.plot(t_sim_QPSS, transient_Vc, "ob",color='green' , label=f'Chutes iniciais Vc')
    plt.plot(t_sim_QPSS, transient_ic1, "ob",color='orange' , label=f'Chutes iniciais IC1')
    plt.plot(t_sim_QPSS, transient_ic2, "ob",color='red' , label=f'Chutes iniciais IC2')
    plt.plot(t_sim_QPSS, results_vc, "ob", label=f'H = {h}')
    plt.plot(t_plot, results_vb)
    plt.title('Tensão no capacitor C1')
    plt.ylabel('(V)')
    plt.xlabel('Tempo (segundos)')
    plt.grid()
    plt.legend()
plt.show()
