import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from math import pi

"""
This code have as objective the transient analysis of a non linear circuit.
To calculate the transient analysis, it is necessary to provide the initial
conditions in the components of circuits and also the equations of circuit.
The simulation depends on previous conditions in each new interation
"""

# inserir as CIs
Ra = 10 ** 3
RL = 50
Vc10 = -1.96951566e+00
Vc20 = -1.89360449e-06
# fonte do Ex1
f1 = 1 * 10 ** 9
f2 = 1.1 * 10 ** 9
# deltat = 1 / (100 * f)
Vm1 = 5
Vm2 = 0.2
C1 = 10 ** (-11)
C2 = 10 ** (-6)


# The matrix from the mesh analysis in circuit. Analysis in t0
def func(x):
    return [Ra * x[0] + Vc10,
            -Vc10 + Vc20 + x[1] * RL]


y = fsolve(func, np.zeros(2))

# previous conditions to new interation
y1 = y[0]
y2 = y[1]

result_vc1_transi = []  # vector to storage the voltage at capacitor 1, to SM
result_vc2_transi = []  # vector to storage the voltage at capacitor 1, to SM

result_vc1 = []  # vector to storage the voltage at capacitor 1
result_vc2 = []  # vector to storage the voltage at capacitor 2
result_ic1 = []  # vector to storage the current at capacitor 1
result_ic2 = []  # vector to storage the current at capacitor 2
nonlinear_element = []

# correntes nos capacitores
ic10 = y1 - y2
ic20 = y2
# storage the values
result_vc1_transi.append(Vc10)
result_vc2_transi.append(Vc20)
result_ic1.append(ic10)
result_ic2.append(ic20)

# notal tensions to the new interation
Va = Vm1 * np.sin(2 * pi * f1 * 0) + Vm2 * np.sin(2 * pi * f2 * 0)
Vb = Vc10
Vc = Vb - Vc20
nonlinear_element.append(((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5))))

deltat = 1 / (100 * f1)
t_sim = np.arange(0, 1 * 1/abs(f1) + deltat, deltat)

# n = int(100 * f1 / f2)
# (t_sim, deltat) = np.linspace(0, 10 * (1 / f1), n, retstep=True)
for t in np.delete(t_sim, 0):
    Vs = Vm1 * np.sin(2 * pi * f1 * t) + Vm2 * np.sin(2 * pi * f2 * t)  # voltage in source in actual time

    # system of nodal analysis to solve in actual time
    def func(x):
        return [x[0] - Vs,
                ((0.1 * np.sign(x[0])) / ((1 + (1.8 / abs(x[0])) ** 5) ** (1 / 5))) - x[1] / Ra - (
                        (2 * C1 / deltat) * (x[1] - Vc10) - ic10) - (
                        (2 * C2 / deltat) * ((x[1] - x[2]) - Vc20) - ic20),
                -((2 * C2 / deltat) * ((x[1] - x[2]) - Vc20) - ic20) + x[2] / RL]


    y = fsolve(func, np.array([Va, Vb, Vc]))  # solving the system of nonlinear equations in t

    # nodal voltages
    Va = y[0]
    Vb = y[1]
    Vc = y[2]
    nonlinear_element_calc = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5)))

    # previous values to next time interation
    ic10 = ((2 * C1 / deltat) * (Vb - Vc10) - ic10)
    ic20 = ((2 * C2 / deltat) * ((Vb - Vc) - Vc20) - ic20)
    Vc10 = Vb
    Vc20 = Vb - Vc

    # storage the values
    nonlinear_element.append(nonlinear_element_calc)
    result_vc1_transi.append(Vc10)
    result_vc2_transi.append(Vc20)
    result_ic1.append(ic10)
    result_ic2.append(ic20)

# plt.plot(t_sim, nonlinear_element)
# plt.title('Tensão no Capacitor 1')
# plt.ylabel('(V)')
# plt.xlabel('Tempo (mili segundos)')
# plt.grid()
# plt.show()
###
# plt.plot (t_plot, result_vc2)
# plt.title ('Tensão no Capacitor 2')
# plt.ylabel ('(V)')
# plt.xlabel ('Tempo (mili segundos)')
# plt.grid()
# plt.show ()
#
# plt.plot (t_plot, result_ic1)
# plt.title ('Corrente da fonte controlada')
# plt.ylabel ('(A)')
# plt.xlabel ('Tempo (mili segundos)')
# plt.grid()
# plt.show ()
# plt.plot (t_plot, result_ic2)
# plt.title ('Corrente no Capacitor 2')
# plt.ylabel ('(A)')
# plt.xlabel ('Tempo (mili segundos)')
# plt.grid()
# plt.show ()
