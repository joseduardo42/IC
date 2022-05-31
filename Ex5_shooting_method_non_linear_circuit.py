import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from math import pi
# from Ex3_transien_nonlinear_circuit import t_plot__aux, result_vc1_transi, aux, result_vc2_transi

"""
This code have as objective the shooting method analysis of a non linear
circuit. To calculate the harmonic balance, it is necessary to provide the initial
conditions in the components of circuits and also the equations of circuit.
The transient analysis is made with the initial conditions obtained with the
shooting method, so the transient is the steady-state
"""

# ICs
Ra = 10 ** 3
RL = 50
C1 = 10 * 10 ** (-12)
C2 = 10 ** (-6)
# Ativ1 source
f = 1 * 10 ** 9
A = 5
deltat = 1 / (100 * f)
tf = (1 / f)
final_resnorm = 0

# time domain
t_sim = np.arange(deltat, tf, deltat)
t_plot = np.arange(0, tf, deltat)
result_vc1 = np.zeros(len(t_plot))
result_vc2 = np.zeros(len(t_plot))
result_ic1 = np.zeros(len(t_plot))
result_ic2 = np.zeros(len(t_plot))


# shooting method
def shooting_method_non_linear(z):
    Vc10 = z[0]
    Vc20 = z[1]

    def func(x):
        return [Ra * x[0] + Vc10,
                -Vc10 + Vc20 + x[1] * RL]

    a = fsolve(func, [0, 0])
    # variables for interactions
    y1 = float(a[0])
    y2 = float(a[1])

    ic10 = y1 - y2  # current C1 i(0)
    ic20 = y2  # current C2 i(0)

    # nodal voltages
    Va = A * np.sin(2 * pi * f * 0)
    Vb = Vc10
    Vc = Vb - Vc20
    result_vc1[0] = Vc10
    result_vc2[0] = Vc20
    result_ic1[0] = ic10
    result_ic2[0] = ic20

    # transient
    i = 1
    for t in t_sim:

        Vs = A * np.sin(2 * pi * f * t)

        # non-linear system to solve in each t
        def func(x):

            return [x[0] - Vs,
                    ((0.1 * np.sign(x[0])) / ((1 + (1.8 / abs(x[0])) ** 5) ** (1 / 5))) - x[1] / Ra - (
                                (2 * C1 / deltat) * (x[1] - Vc10) - ic10) - (
                                (2 * C2 / deltat) * ((x[1] - x[2]) - Vc20) - ic20),
                    ((2 * C2 / deltat) * ((x[1] - x[2]) - Vc20) - ic20) - x[2] / RL]

        y = fsolve(func, [Va, Vb, Vc], full_output=True)


        # initial guess for the next interaction
        Va = float(y[0][0])
        Vb = float(y[0][1])
        Vc = float(y[0][2])

        # valtage and current in each capacitor for the next iteraction
        ic10 = ((2 * C1 / deltat) * (Vb - Vc10) - ic10)
        ic20 = ((2 * C2 / deltat) * ((Vb - Vc) - Vc20) - ic20)
        Vc10 = Vb
        Vc20 = Vb - Vc

        # capacitors voltage
        result_vc1[i] = Vc10
        result_vc2[i] = Vc20
        result_ic1[i] = ic10
        result_ic2[i] = ic20
        i += 1

    # discretized equations
    return [result_vc1[-1] - z[0],
            result_vc2[-1] - z[1]]


# solving shooting method
sm_nonlin_result = fsolve(shooting_method_non_linear, [0, 0])
print(sm_nonlin_result)
# plotting 1 cycle of transient in C1
plt.plot(t_plot, result_vc1)
#plt.plot(t_plot__aux, result_vc1_transi[-aux:])
plt.title('Tensão no capacitor 1')
plt.ylabel('Tensão no capacitor 1 (V)')
plt.xlabel('Tempo (segundos)')
plt.grid()
plt.show()
#
# plotting 1 cycle of transient in C2
# plt.plot (t_plot, result_vc2)
# plt.plot (t_plot__aux, result_vc2_transi[-aux:])
# plt.title ('Tensão no capacitor 2')
# plt.ylabel ('Tensão no capacitor (V)')
# plt.xlabel ('Tempo (segundos)')
# plt.grid ()
# plt.show ()
