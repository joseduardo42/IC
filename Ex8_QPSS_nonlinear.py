import numpy as np
from numpy import linalg, pi
import matplotlib.pyplot as plt
from math import sin, cos, pi
from scipy.optimize.minpack import fsolve
from Ex7_superposition import t_sim, results_vb
# ICs
C1 = 10 * 10 ** -12
C2 = 1 * 10 ** -6
R1 = 1 * 10 ** 3
RL = 50
f1 = 1 * 10 ** 9
f2 = 1.01 * 10 ** 9
fA = (f1 + f2)/2
fB = (abs(f1 - f2))
w1 = 2 * pi * f1 
w2 = 2 * pi * f2
Vm1 = 5
Vm2 = 0.2
h1 = 5

T1 = (1 / f1)
T = (1 / f2)
T1_nl = (1 / fA)
T_nl = (1 / fB)

for h in range(1, 2):
    k = 2 * h + 1
    k1 = 2 * h1 + 1
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
    shooting_Vb = np.zeros(k)
    shooting_Vc = np.zeros(k)

    transient_Vb = np.zeros(k)
    transient_Vc = np.zeros(k)    

    def qpss(shooting_voltage):
        
        for i in range(k):
            
            N = k1
            (t_sim, deltat) = np.linspace(i * (T_nl / k), i * (T_nl / k) + T1_nl, N, retstep=True)            
            
            Va = Vm1 * np.sin(2 * pi * f1 * t_sim[0]) + Vm2 * np.sin(2 * pi * f2 * t_sim[0])
    
            Vc10 = shooting_voltage[i]
            Vc20 = shooting_voltage[i] - shooting_voltage[k + i]

            # system of nodal analysis to solve in actual time
            def func(x):

                if Va == 0:
                    i_fnl = 0
                else:
                    i_fnl = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5)))

                return [i_fnl * R1 - x[0] * R1 + Vc10,
                        -Vc10 + Vc20 + x[1] * RL]


            initial_conditions = fsolve(func, np.array(
                        [0, 0]))  # solving the system of nonlinear equations in t

            ic10 = float(initial_conditions[0]) - float(initial_conditions[1])
            ic20 = float(initial_conditions[1])

            for t in np.delete(t_sim, 0):
                Vs = Vm1 * np.sin(2 * pi * f1 * t) + Vm2 * np.sin(2 * pi * f2 * t)  # voltage in source in actual time

                # system of nodal analysis to solve in actual time
                def func(x):

                    if  Vs == 0:
                        i_fnl = 0
                    else:
                        i_fnl = ((0.1 * np.sign(Vs)) / ((1 + (1.8 / abs(Vs)) ** 5) ** (1 / 5)))

                        return [i_fnl * R1 - x[0] * R1 + ((deltat/(2*C1))*((x[0]-x[1]) + ic10) + Vc10),
                                
                                -((deltat/(2 * C1)) * ((x[0]-x[1]) + ic10) + Vc10) + ((deltat/(2*C2)) * (x[1] + ic20) + Vc20) + x[1] * RL]


                y = fsolve(func, np.array([0, 0]))  # solving the system of nonlinear equations in t

                # nodal voltages
                ic1 = y[0]
                ic2 = y[1]
                Vc10 = (deltat/(2*C1))*((ic1-ic2) + ic10) + Vc10
                Vc20 = (deltat/(2*C2))*(ic2 + ic20) + Vc20
                Vb = Vc10
                Vc = Vc10 - Vc20
                nonlinear_element_calc = ((0.1 * np.sign(Vs)) / ((1 + (1.8 / abs(Vs)) ** 5) ** (1 / 5)))

                # previous values to next time interation
                ic10 = ic1 - ic2
                ic20 = ic2
            
                # transient results
                
                transient_Vb[i] = Vb
                transient_Vc[i] = Vc
        
        return np.concatenate([
            
            transient_Vb - D @ shooting_voltage[0: k] ,
            transient_Vc - D @ shooting_voltage[k:2 * k],

        ])


    # solve qpss function
    amplitudes_guess = np.ones(2*k)
    y = fsolve(qpss, amplitudes_guess)


    t_sim_QPSS = np.array([i * (T_nl / k) for i in range(k)])

    # plot results
    plt.plot(t_sim, results_vb, label=f'PAC')
    plt.plot(t_sim_QPSS,y[0: k], "ob", color="red", label=f'QPSS')
    plt.title('Tensão em Vb, linearização em torno de fB (omega2 no PAC)')
    plt.ylabel('(V)')
    plt.xlabel('Tempo (segundos)')
    plt.grid()
    plt.legend()
plt.show()
