import numpy as np
from numpy.lib.function_base import append
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from math import pi

#inserir as CIs
Ra = 10**3
RL = 50
Vc10 = 10
Vc20 = 20
#fonte do Ex2
Vm1 = 60
Vm2 = 40
f1 = 100
f2 = 8
Vs = 0
deltat = 1/(20 * f1)
tf = 2*(1/4)

C1 = 10**(-11)
C2 = 10**(-6)

#usando fsolve para resolver a primeira analise em t=0
def func(x):
    return[ Ra*x[0] + Vc10,
            -Vc10 + Vc20 + x[1]*RL]

y = fsolve(func, [0, 0])

result_vc1 = [] #vetor que armazena a tensão no capacitor 1
result_vc2 = [] #vetor que armazena a tensão no capacitor 2
result_ic1 = [] #vetor que armazena a corrente no capacitor 1
result_ic2 = [] #vetor que armazena a corrente no capacitor 2

#correntes nos capacitores
ic10 = float (y[0]-y[1])
ic20 = float (y[1])
#guardando valores em um vetor
result_vc1.append (Vc10)
result_vc2.append (Vc20)
result_ic1.append (ic10)
result_ic2.append (ic20)
#variavei que irão ser usadas no for, assumirão outros valores a cada iteração
y1 = ic10
y2 = ic20

t_sim = np.arange(deltat, tf+deltat, deltat)#tempo de simulação
t_plot = np.arange(0, tf+deltat, deltat) #vetor com todas posições a ser plotado

for t in t_sim:

    Vs = (Vm1*np.sin(2*pi*f1*t) + Vm2*np.sin(2*pi*f2*t)) #tensão da fonte no tempo desejado
    #equações da segunda analise
    def func(x):
        return[Ra*(x[0]-((0.1*np.sign(Vs)) / ((1 + (1.8/abs(Vs))**5)**(1/5)))) + Vc10 + (1/(2*C1))*deltat*((x[0] - x[1]) + ic10),
               -(Vc10 + 1/(2*C1)*deltat*((x[0] - x[1]) + ic10)) + (Vc20 + 1/(2*C2)*deltat*(x[1] + ic20)) + RL*x[1]]
    
    y = fsolve(func, [y1, y2])
    #valores das tensões para a próxima iteração
    Vc10 = (Vc10 + 1/(2*C1)*deltat*((y[0] - y[1]) + ic10))
    Vc20 = (Vc20 + 1/(2*C2)*deltat*(y[1] + ic20))
    #atualizando variáveis do chute inicial
    y1 = float (y[0])
    y2 = float (y[1])
    #valores das correntes para a próxima iteração
    ic10 = y1 - y2
    ic20 = y2
    #armazenando variáveis
    result_vc1.append (Vc10)
    result_vc2.append (Vc20)
    result_ic1.append (ic10)
    result_ic2.append (ic20)


plt.plot (t_plot, result_vc1)
plt.title ('Tensão no Capacitor 1')
plt.ylabel ('(V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()

plt.plot (t_plot, result_vc2)
plt.title ('Tensão no Capacitor 2')
plt.ylabel ('(V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()

plt.plot (t_plot, result_ic1)
plt.title ('Corrente no Capacitor 1')
plt.ylabel ('(A)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()

plt.plot (t_plot, result_ic2)
plt.title ('Corrente no Capacitor 2')
plt.ylabel ('(A)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()