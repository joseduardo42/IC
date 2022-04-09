import numpy as np
from numpy import pi, cos, sin
import matplotlib.pyplot as plt
from Ex7_PAC_HB_linear import x1_lin, x2_lin, x3_lin, X_c1_lin
from Ex3_transien_nonlinear_circuit import nonlinear_element as nonlinear_element_transient
from Ex7_PAC_HB_nonlinear import x1_1tom, x2_1tom, x3_1tom, X_c1_1tom, h1

# params
C1 = 10 * 10 ** -12
C2 = 1 * 10 ** -6
R1 = 1 * 10 ** 3
RL = 50
f1 = 1 * 10 ** 9
w1 = 2 * pi * f1
f2 = 1.1 * 10 ** 9
w2 = 2 * pi * f2
Vm1 = 5
Vm2 = 0.2
h2 = int(h1 / 2)
k2 = 2 * (2 * h2 + 1)
T1 = (1 / f1)
T = (1 / f2)

results_va = []
results_vb = []
results_vc = []
nonlinear_element = []

deltat = 1 / (100 * f1)
t_sim = np.arange(0, 10 * 1/f1 + deltat, deltat)
for t in t_sim:
    sinandcos_lin = np.array([f((w2 + j * w1) * t) for j in range(-h2, h2 + 1, 1) for f in (sin, cos)])
    sinandcos_nonlin = np.array([1] + [f((w1 * (j + 1)) * t) for j in range(h1) for f in (sin, cos)])
    Va_time = sinandcos_nonlin @ x1_1tom + sinandcos_lin @ x1_lin
    Vb_time = sinandcos_nonlin @ x2_1tom + sinandcos_lin @ x2_lin
    Vc_time = sinandcos_nonlin @ x3_1tom + sinandcos_lin @ x3_lin
    Vc1_time = sinandcos_nonlin @ X_c1_1tom + sinandcos_lin @ X_c1_lin

    dependent_source = Vc_time / RL + Vb_time / R1 + Vc1_time

    results_va.append(Va_time)
    results_vb.append(Vb_time)
    results_vc.append(Vc_time)
    nonlinear_element.append(dependent_source)

# plot results
plt.plot(t_sim, nonlinear_element, label='HB')
plt.legend(loc="upper right")
plt.plot(t_sim, nonlinear_element_transient, label='Transitório')
plt.legend(loc="upper right")
plt.title('Corrente da fonte controlada')
plt.ylabel('(A)')
plt.xlabel('Tempo (mili segundos)')
plt.grid()
plt.show()
