import sympy as sym
import numpy as np
import math
import control as ctrl
import matplotlib.pyplot as plt

m, g, d, δ, r, R, L, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i, y, x1, Xe, x2 = sym.symbols('m, g, d, δ, r, R, L, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i, y, x1, Xe, x2')

F = (v**2)*(2*c/m*3*(R*y+L)**2) + (2/3)*g*sym.sin(φ) - (2/3)*b*x2 - (2/3*m)*k1*(x1 - d)-(2/3*m)*k2*(x1-d)**3

Xmax = δ
Xmin = d + (m * g * sym.sin(φ)/k)

# Equilibrium point
Feq = 0
x2eq = 0
x1eq = Xe = (0.25 * Xmax) + (0.75 * Xmin)
#veq = ieq * R

#dphi_F_eq = dphi_F.subs([(F, Feq), (x2, x2eq), (x4, x4eq)])
#dphi_x3_eq = dphi_x3.subs([(F, Feq), (x3, x3eq), (x4, x4eq)])

q = 2 * c / 3 * m * (R * y + L) ** 2
u = 2 * b / 3 * m
w = 2 * (k1 + k2 * 3 * (x1eq - d) ** 2) / 3 * m


#Fmag = c(i**2/y**2)

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

def evaluate_at_given_parameters(z):
    """
    :param z:
    :return:
    """
    return float(z.subs( [(m, m_value), (g, g_value), (d, d_value), (δ, δ_value), (r, r_value), (R, R_value), (L0, L0_value), (L1, L1_value), (α, α_value), (c, c_value), (k, k_value), (b, b_value), (φ, φ_value)]))

L= L0 + (L1 * math.exp(-α_value * (α_value - x1)))

#Fmag_value = evaluate_at_given_parameters(Fmag)
#L_value = evaluate_at_given_parameters(L)
#y_value = evaluate_at_given_parameters(y)

#F_value = evaluate_at_given_parameters(F)

q, u, w = sym.symbols('q, u, w', real=True, positive=True)
s, t = sym.symbols('s, t')

transfer_function = ctrl.TransferFunction([q], [1, u, w])

def pid (kp, ki, kd):

    diff    = ctrl.TransferFunction ([1, 0], 1)
    intgr   = ctrl.TransferFunction (1, [1, 0])
    pid_f = kp + kd * diff + ki * intgr

    return pid_f

kp =
ki =
kd =

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

