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
Vc10 = -73.29215789
Vc20 = -76.95672331
# fonte do Ex1
f = 10 ** 9
deltat = 1 / (100 * f)
tf = (1 / f)
A = 100

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

# correntes nos capacitores
ic10 = y1 - y2
ic20 = y2
# storage the values
result_vc1_transi.append(Vc10)
result_vc2_transi.append(Vc20)
result_ic1.append(ic10)
result_ic2.append(ic20)

# notal tensions to the new interation
Va = A * np.sin(2 * pi * f * 0)
Vb = Vc10
Vc = Vb - Vc20

t_sim = np.arange(deltat, tf, deltat)  # time vector to simulation, without t0
t_plot = np.arange(0, tf, deltat)  # vector to plot in each time of simulation

for t in t_sim:
    Vs = A * np.sin(2 * pi * f * t)  # voltage in source in actual time

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

    # previous values to next time interation
    ic10 = ((2 * C1 / deltat) * (Vb - Vc10) - ic10)
    ic20 = ((2 * C2 / deltat) * ((Vb - Vc) - Vc20) - ic20)
    Vc10 = Vb
    Vc20 = Vb - Vc

    # storage the values
    result_vc1_transi.append(Vc10)
    result_vc2_transi.append(Vc20)
    result_ic1.append(ic10)
    result_ic2.append(ic20)

# variables for plot comparing to shooting method. Otherwise, use t_plot and result_vc
t_plot__aux = np.arange(0, 1 / f, deltat)
aux = len(t_plot__aux)

plt.plot(t_plot__aux, result_vc1_transi[-aux:])
plt.title('Tensão no Capacitor 1')
plt.ylabel('(V)')
plt.xlabel('Tempo (mili segundos)')
plt.grid()
plt.show()
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
