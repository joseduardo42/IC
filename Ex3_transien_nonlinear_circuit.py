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
Vc10 = -1.452485044906209
Vc20 = -1.4350046838156988e-05
# fonte do Ex1
f1 = 1.0 * 10 ** 9
f2 = 1.01 * 10 ** 9
# deltat = 1 / (100 * f)
Vm1 = 1
Vm2 = 1
C1 = 10 * 10 ** (-12)
C2 = 1 *  10 ** (-6)

T1_nl = (1/(abs(f1 + f2)/2))
k1 = 2 * 10000 + 1

N = int(k1)
(t_sim, deltat) = np.linspace(0, (1/abs(f1 - f2)), N, retstep=True)


result_vc1_transi = []  # vector to storage the voltage at capacitor 1, to SM
result_vc2_transi = []  # vector to storage the voltage at capacitor 1, to SM

result_vc1 = []  # vector to storage the voltage at capacitor 1
result_vc2 = []  # vector to storage the voltage at capacitor 2
result_ic1 = []  # vector to storage the current at capacitor 1
result_ic2 = []  # vector to storage the current at capacitor 2
nonlinear_element = []

# notal tensions to the new interation
Va = Vm1 * np.sin(2 * pi * f1 * t_sim[0]) + Vm2 * np.sin(2 * pi * f2 * t_sim[0] )
Vb = Vc10
Vc = Vb - Vc20

if  Va == 0:
    nonlinear_element_calc = 0
else:
    nonlinear_element_calc = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5)))

def func(x):

    if  Va == 0:
        i_fnl = 0
    else:
        i_fnl = ((0.1 * np.sign(Va)) / ((1 + (1.8 / abs(Va)) ** 5) ** (1 / 5)))

    return [i_fnl * Ra - x[0] * Ra + Vc10,
            -Vc10 + Vc20 + x[1] * RL]

initial_conditions = fsolve(func, np.array(
                [0, 0]))  # solving the system of nonlinear equations in t

# correntes nos capacitores
ic10 = float(initial_conditions[0]) - float(initial_conditions[1])

ic20 = float(initial_conditions[1])
# storage the values
result_vc1_transi.append(Vb)
result_vc2_transi.append(Vc20)
result_ic1.append(ic10)
result_ic2.append(ic20)

for t in np.delete(t_sim, 0):
    
    Vs = Vm1 * np.sin(2 * pi * f1 * t) + Vm2 * np.sin(2 * pi * f2 * t) # voltage in source in actual time

    # system of nodal analysis to solve in actual time
    def func(x):

        if  Vs == 0:
            i_fnl = 0
        else:
            i_fnl = ((0.1 * np.sign(Vs)) / ((1 + (1.8 / abs(Vs)) ** 5) ** (1 / 5)))

            return [i_fnl * Ra - x[0] * Ra + ((deltat/(2*C1))*((x[0]-x[1]) + ic10) + Vc10),
                    
                    -((deltat/(2 * C1)) * ((x[0]-x[1]) + ic10) + Vc10) + ((deltat/(2*C2)) * (x[1] + ic20) + Vc20) + x[1] * RL]


    y = fsolve(func, np.array([0, 0]))  # solving the system of nonlinear equations in t

    # nodal voltages
    ic1 = y[0]
    ic2 = y[1]
    Vc10 = (deltat/(2*C1))*((ic1-ic2) + ic10) + Vc10
    Vc20 = (deltat/(2*C2))*(ic2 + ic20) + Vc20

    if  Vs == 0:
        nonlinear_element_calc = 0
    else:
        nonlinear_element_calc = ((0.1 * np.sign(Vs)) / ((1 + (1.8 / abs(Vs)) ** 5) ** (1 / 5)))
    # previous values to next time interation
    ic10 = ic1 - ic2
    ic20 = ic2
    

    # storage the values
    nonlinear_element.append(nonlinear_element_calc)
    result_vc1_transi.append(Vc10)
    result_vc2_transi.append(Vc20)
    result_ic1.append(ic10)
    result_ic2.append(ic20)

print(Vc10, Vc20)
# plt.plot(t_sim, result_vc1_transi)
# plt.title('Tensão no Capacitor 1')
#plt.ylabel('(V)')
#plt.xlabel('Tempo (mili segundos)')
#plt.grid()
#plt.show()
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
