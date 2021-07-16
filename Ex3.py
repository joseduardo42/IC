import numpy as np
from numpy.lib.function_base import append
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#inserir as CIs
Ra = 10**3
RL = 50
Vc10 = 10
Vc20 = 20
Vs = 0
w = 100000
A = 1.2
deltat = 1/(20 * w)
tf = 5*1/w
C1 = 10**(-11)
C2 = 10**(-6)

#construindo matriz a partir da análise de malhas
def func(x):
    return[ Ra*x[0] + Vc10,
            -Vc10 + Vc20 + x[1]*RL]

y = fsolve(func, [0, 0])

result_vc1 = [] #vetor que armazena a tensão no capacitor
result_vc2 = [] #vetor que armazena a tensão no capacitor
result_ic1 = [] #vetor que armazena a corrente na fonte
result_ic2 = [] #vetor que armazena a corrente na fonte

ic10 = float (y[0]-y[1])
ic20 = float (y[1])

result_vc1.append (Vc10)
result_vc2.append (Vc20)
result_ic1.append (ic10)
result_ic2.append (ic20)

y1 = 0
y2 = ic10
y3 = ic20

t_sim = np.arange(deltat, tf+deltat, deltat)
t_plot = np.arange(0, tf+deltat, deltat) #vetor com todas posições a ser plotado

for t in t_sim:

    Vs = (100*np.sin(w*t)) #tensão da fonte no tempo desejado

    def func(x):
        return[((0.1*np.sign(Vs)) / (1 + (1.8/abs(Vs))**5)**(1/5)) + Ra*(x[0]-x[1]),
                Ra*(x[1]-x[0]) + Vc10 + (1/(2*C1))*deltat*((x[1] - x[2]) + ic10),
                -(1/(2*C1))*deltat*((x[1] - x[2]) + ic10) + (1/(2*C2))*deltat*(x[2] + ic20)]
    
    y = fsolve(func, [y1, y2, y3])

    Vc10 = (1/(2*C1))*deltat*((y[1] - y[2]) + ic10)
    Vc20 = (1/(2*C2))*deltat*(y[2] + ic20)
    
    y1 = float (y[0])
    y2 = float (y[1])
    y3 = float (y[2])

    ic10 = y2 - y3
    ic20 = y3

    result_vc1.append (Vc10)
    result_vc2.append (Vc20)
    result_ic1.append (ic10)
    result_ic2.append (ic20)

'''
plt.plot (t_plot, result_ic2)
plt.title ('Tensão no capacitor')
plt.ylabel ('Tensão no capacitor (V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()
'''