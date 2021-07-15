#import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg, pi

#inserir as CIs
Ra = 1000
RL = 50
Vc10 = 10
Vc20 = 20
Vs = 0
#V1 = 0; V2 = 0; V3 = 0


#construindo matriz a partir da análise de malhas
A1 = np.array([[1, 0, 0],
              [0, 1/Ra, 0],
              [0, 0, 1/RL]], np.float64)
              
b = np.array([0, -(Vc10+Vc20), Vc20], np.float64)

x = linalg.solve(A1, b)

print (x)