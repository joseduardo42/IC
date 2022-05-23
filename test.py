from turtle import color
import numpy as np
from numpy import linalg, pi
import matplotlib.pyplot as plt
from math import sin, cos, pi

from scipy.optimize.minpack import fsolve

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
h = 2
final_resnorm = 0

T1 = (1 / f1)
T = (1 / f2)
T1_nl = (1 / fA)
T_nl = (1 / fB)


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
shooting_Vb = np.zeros(k)
shooting_Vc = np.zeros(k)

transient_Va = np.zeros(k)
transient_Vb = np.zeros(k)
transient_Vc = np.zeros(k)



def qpss(shooting_voltage):

    # vector of unknowns

    for i in range(k):
        shooting_Vb[i] = shooting_voltage[i]
        shooting_Vc[i] = shooting_voltage[k + i]


    ################## shooting #####################
    ################## transient (5) #####################

    for i in range(k):

        deltat = T1_nl / k
        t_sim = np.arange(i * (T_nl / k), i * (T_nl / k) + T1_nl, deltat)
        
        
        Va = Vm1 * np.sin(2 * pi * f1 * t_sim[0]) + Vm2 * np.sin(2 * pi * f2 * t_sim[0])
        
        Vc10 = shooting_Vb[i]
        Vc20 = shooting_Vb[i] - shooting_Vc[i]
        
        # system of nodal analysis to solve in actual time
        def func(x):

            if Va == 0:
                i_fnl = 0
            else:
                i_fnl = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5)))

            return [i_fnl * R1 - x[0] * R1 + Vc10,
                    -Vc10 + Vc20 + x[1] * RL]


        initial_conditions = fsolve(func, [0, 0])  # solving the system of nonlinear equations in t

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
            
            transient_Vb[i] = Vb
            transient_Vc[i] = Vc

    
    return np.concatenate([
        
        transient_Vb - D @ shooting_Vb,
        transient_Vc - D @ shooting_Vc,

    ])


# solve qpss function
amplitudes_guess = np.zeros(k * 2)

y = fsolve(qpss, amplitudes_guess)

QPSS_result_Vb = y[0: k]
QPSS_result_Vc = y[k: 2 * k]

results_transient = []
second_transient_Vb = np.zeros(k)

for i in range(k):
    deltat = T1_nl / k
    t_sim = np.arange(i * (T_nl / k), i * (T_nl / k) + T1_nl, deltat)
    
    Va = Vm1 * np.sin(2 * pi * f1 * t_sim[0]) + Vm2 * np.sin(2 * pi * f2 * t_sim[0])

    Vc10 = QPSS_result_Vb[i]
    Vc20 = QPSS_result_Vb[i] - QPSS_result_Vc[i]
    
    results_transient.append(QPSS_result_Vb[i])
    
    # system of nodal analysis to solve in actual time
    def func(x):

        if Va == 0:
            i_fnl = 0
        else:
            i_fnl = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5)))

        return [i_fnl * R1 - x[0] * R1 + Vc10,
                -Vc10 + Vc20 + x[1] * RL]


    initial_conditions = fsolve(func, [0, 0])  # solving the system of nonlinear equations in t

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
        
        results_transient.append(Vb)

        # capacitor current
        ic10 = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5))) - Vb / R1 - Vc / RL
        ic20 = Vc / RL

        # capacitor voltage
        Vc10 = Vb
        Vc20 = Vb - Vc
        
        second_transient_Vb[i] = Vb
        
transient_1 = results_transient[0: k]
transient_2 = results_transient[k: 2*k]
transient_3 = results_transient[2*k: 3*k]
transient_4 = results_transient[3*k: 4*k]
transient_5 = results_transient[4*k: 5*k]

second_transient_Vb_delay = D @ second_transient_Vb

plt.plot(t_sim, transient_1, "ob",color='red' , label=f'Primeiro transitório')
plt.plot(t_sim, transient_2, "ob",color='yellow' , label=f'Segundo transitório')
plt.plot(t_sim, transient_3, "ob",color='green' , label=f'Terceiro transitório')
plt.plot(t_sim, transient_4, "ob",color='orange' , label=f'Quarto transitório')
plt.plot(t_sim, transient_5, "ob",color='brown' , label=f'Quinto transitório')
plt.plot(t_sim, second_transient_Vb_delay, "ob",color='pink' , label=f'Ultimo ponto de cada transitorio')
plt.plot(t_sim, second_transient_Vb, "ob",color='blue' , label=f'D @ ultimo ponto de cada transitorio')
plt.title('Tensão no capacitor C1')
plt.ylabel('(V)')
plt.xlabel('Tempo (segundos)')
plt.grid()
plt.legend()
plt.show()
