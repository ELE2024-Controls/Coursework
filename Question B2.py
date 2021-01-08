import sympy as sym
import numpy as np
import math
import control as ctrl
import matplotlib.pyplot as plt

m, g, d, δ, r, R, L, L0, L1, α, c, k, b, φ, veq, xe, x, x1, x2, y, v = sym.symbols('m, g, d, δ, r, R, L, L0, L1, α, c, k, b, φ, veq, xe, x, x1, x2, y, v')

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

# defining x1 x2 and dx2

x1 = x
x2 = sym.diff(x1)
dx2 = sym.diff(x2)

# Equilibrium point
dx2_eq = 0
x2eq = 0

X_min = d_value + (m_value * g_value * sym.sin(φ_value)/k_value)
X_max = δ_value
xe = 0.75 * X_min + 0.25 * X_max

y_value = δ_value - xe
L_value = L0_value + (L1_value * sym.exp(-α_value * y_value))


dx2_eq = 2*m_value*c_value*(veq**2)/(R_value*y_value+L_value)**2 + 2*m_value**2*g_value*sym.sin(φ_value) - 2*m_value*k_value*xe - 2*m_value*b_value*x2eq

print("")
print("xe value:")
print("")
sym.pprint(xe)
print("")

#-----------code to linearise system----------------




#---------------------------------------------------

print("")
print("linearised equation:")
print("")
#sym.pprint(linearised_eqn)
print("")
