import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg, pi

#inserir as CIs
R1 = 1
R2 = 2
R3 = 3
C = 1
Vm1 = 60
Vm2 = 40
f1 = 100
f2 = 8
deltat = 1/(20 * f1)
tf = 2*(1/4)
vs = 0
vc0 = 20

#construindo matriz a partir da análise de malhas
A1 = np.array([[R1, -R1],
            [-R1, (R1 + R2 + R3)]], np.float64)
b = np.array([[vs], [-vc0]], np.float64)

x = linalg.solve(A1,b)

i0 = float (x[1]) #corrente a ser utilizada no próximo equacionamento, instante anterior

t_sim = np.arange(deltat, tf+deltat, deltat) #vetor com cada tempo a ser utilizado. Iniciei em t0+deltat, pois t0 já foi feito a análise anteriormente
t_plot = np.arange(0, tf+deltat, deltat) #vetor com todas posições a ser plotado
result_vc = [] #vetor que armazena a tensão no capacitor
result_ic = [] #vetor que armazena a corrente na fonte
result_is = [] #vetor que armazena a corrente da fonte
#armazenar os valores de t0
result_vc.append (vc0)
result_ic.append (i0)
result_is.append (float (x[0]))

for t in t_sim:

        vs = (Vm1*np.sin(2*pi*f1*t) + Vm2*np.sin(2*pi*f2*t)) #tensão da fonte no tempo desejado

        #sistema a ser resolvido a cada instante t        
        A2 = np.array([[R1, -R1, 0],
            [-R1, (R1 + R2 + R3), 1],
            [0, -deltat/(2*C), 1]], np.float64)
        b = np.array([[vs], [0], [vc0 + i0*deltat/(2*C)]], np.float64)

        x = linalg.solve(A2, b) #sistema algébrico linear a ser calculado a cada t
        
        #armazenar corrente calculada da fonte
        i_s = float (x[0])
        result_is.append (i_s)
        #corrente no capacitor a ser usada na próxima iteração
        ic = float (x[1])
        result_ic.append (ic)
        i0 = ic
        #tensão no capacitor a ser usada na próxima iteração
        vc = float (x[2])
        result_vc.append (vc)
        vc0 = vc
        
        

plt.plot (t_plot, result_vc)
plt.title ('Tensão no capacitor')
plt.ylabel ('Tensão no capacitor (V)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()

plt.plot (t_plot, result_ic)
plt.title ('Corrente no capacitor')
plt.ylabel ('Corrente no capacitor (A)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()

plt.plot (t_plot, result_is)
plt.title ('Corrente da fonte')
plt.ylabel ('Corrente da fonte (A)')
plt.xlabel ('Tempo (milisegundos)')
plt.grid()
plt.show ()
