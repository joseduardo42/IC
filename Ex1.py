import numpy as np
from numpy import linalg

#inserir as CIs
R1 = 1
R2 = 2
R3 = 3
C = 1
A = 100
w = 100
deltat = 1/(20 * w)
tf = 5*1/w
vs = 0
vc0 = 20

#construindo matriz a partir da análise de malhas
A1 = np.array([[R1, -R1],
            [-R1, (R1 + R2 + R3)]], np.float64)
b = np.array([[vs], [-vc0]], np.float64)

x = linalg.solve(A1,b)

i0 = float (x[1]) #corrente a ser utilizada no próximo equacionamento, instante anterior

t_sim = np.arange(deltat, tf+deltat, deltat) #vetor com cada tempo a ser utilizado. Iniciei em t0+deltat, pois t0 já foi feito a análise anteriormente
m = t_sim.shape[0]
result_vc = np.zeros(m+1) #vetor que armazena a tensão no capacitor
result_ic = np.zeros(m+1) #vetor que armazena a corrente na fonte

#armazenar os valores de t0
result_vc[0] = vc0
result_ic[0] = i0
i=0
for t in t_sim:
        #contador para definir a posição dos itens no vetor
        i += 1
        vs = (100*np.sin(w*t)) #tensão da fonte no tempo desejado

        #sistema a ser resolvido a cada instante t        
        A2 = np.array([[R1, -R1, 0],
            [-R1, (R1 + R2 + R3), 1],
            [0, -deltat/(2*C), 1]], np.float64)
        b = np.array([[vs], [0], [vc0 + i0*deltat/(2*C)]], np.float64)

        x = linalg.solve(A2, b) #sistema algébrico linear a ser calculado a cada t
        
        #corrente no capacitor a ser usada na próxima iteração
        ic = float (x[1])
        result_ic[i] = ic
        i0 = ic
        #tensão no capacitor a ser usada na próxima iteração
        vc = float (x[2])
        result_vc[i] = vc
        vc0 = vc

print (result_vc)
print (result_ic)
