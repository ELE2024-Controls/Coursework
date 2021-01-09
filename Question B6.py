import sympy as sym
import numpy as np
import math
import control as ctrl
import matplotlib.pyplot as plt

m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, x, xe, x1, x2, k1, k2, ieq, y, veq = sym.symbols('m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, x, xe, x1, x2, k1, k2, ieq, y, veq')

# GIVEN VALUES!
m_value = 0.425
g_value = 9.81
d_value = 0.42
δ_value = 0.65
r_value = 0.125
R_value = 53
L0_value = 0.120
L1_value = 0.025
α_value = 1.2
c_value = 6815
k_value = 1880
b_value = 10.4
φ_value = 42

x1 = xe
x2 = sym.diff(x1)
dx2 = sym.diff(x2)

y = δ - xe
L = L0 + (L1 * sym.exp(-α * y))

dx2 = 2*c*v**2/(R*y+L)**2 + 2*g*sym.sin(φ) - 2*k*x1/m - 2*b*x2/m

# defining Xe(x1) and X2 at Equilibrium point

x2eq = 0
dx2eq = 0

X_min = d_value + (m_value * g_value * sym.sin(φ_value)/k_value)
X_max = δ_value
xe = 0.75 * X_min + 0.25 * X_max

y_value = δ_value - xe
L_value = L0_value + (L1_value * sym.exp(-α_value * y_value))


A, B, C = sym.symbols('A, B, C', real=True, positive=True)
s, t = sym.symbols('s, t')

A = 2*c*2*v/m*(R*y*L)**2

B = 2*k/m

C = 2*b/m

transfer_function = ctrl.TransferFunction([A], [1, C, B])

def pid (kp, ki, kd):

    diff    = ctrl.TransferFunction ([1, 0], 1)
    intgr   = ctrl.TransferFunction (1, [1, 0])
    pid_f = kp + kd * diff + ki * intgr

    return pid_f

kp = 1
ki = 1
kd = 1

controller = pid(kp, ki, kd)

# 1 sample every 30ms
sampling_rate = 33.3333
n_points = 500

closed_loop_f = ctrl.feedback(transfer_function, controller)
t, y = ctrl.impulse_response(closed_loop_f, sampling_rate)

plt.plot(t, y)
plt.suptitle('System PID controller', fontsize=14)
plt.xlabel('time(s)')
plt.ylabel('distance from Xsp (m)')
plt.grid()
plt.show()

