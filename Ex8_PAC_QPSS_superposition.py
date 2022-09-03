import numpy as np
from numpy import pi, cos, sin
import matplotlib.pyplot as plt
from Ex8_PAC_QPSS_linear import x1_lin, x2_lin, x3_lin, X_c1_lin, f2
from Ex8_PAC_QPSS_nonlinear import f1, h1
from Ex8_QPSS_CONFER import t_sim_QPSS, x2_1tom
from Ex3_transien_nonlinear_circuit import result_vc1_transi, t_sim

# params
C1 = 10 * 10 ** -12
C2 = 1 * 10 ** -6
R1 = 1 * 10 ** 3
RL = 50
w1 = 2 * pi * f1
w2 = 2 * pi * f2
h2 = 1
k2 = 2 * (2 * h2 + 1)

results_va = []
results_vb = []
results_vc = []
nonlinear_element = []

for t in t_sim_QPSS:
    sinandcos_lin = np.array([f((w2 + j * w1) * t) for j in range(-h2, h2 + 1, 1) for f in (sin, cos)])
    Va_time = sinandcos_lin @ x1_lin
    Vb_time = sinandcos_lin @ x2_lin
    Vc_time = sinandcos_lin @ x3_lin
    Vc1_time = sinandcos_lin @ X_c1_lin
    dependent_source = Vc_time / RL + Vb_time / R1 + Vc1_time

    results_va.append(Va_time)
    results_vb.append(Vb_time)
    results_vc.append(Vc_time)
    nonlinear_element.append(dependent_source)

# plot results
plt.plot(t_sim, result_vc1_transi, label=f'Transitorio')
plt.plot(t_sim_QPSS, results_vb, "o", label='PAC')
plt.plot(t_sim_QPSS, x2_1tom, "o",  label=f'QPSS')
plt.legend(loc="upper right")
plt.title('Tensão em Vb, linearização em torno de omega1 (fA no QPSS)')
plt.ylabel('(V)')
plt.xlabel('Tempo (mili segundos)')
plt.grid()
plt.show()
