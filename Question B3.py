import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, xe, x1, x2, k1, k2, i, y, veq, L = sym.symbols('m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, xe, x1, x2, k1, k2, i, y, veq, L')

# Question B3 Determine and plot the impulse and step responses of the transfer function of the linearised system

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

# finding q, u and w, laying out transfer function

A = 2*c*2*v/m*(R*y*L)**2

print ("")
print ("A:")
print ("")
sym.pprint(A)
print ("")

B = 2*k/m

print ("")
print ("B:")
print ("")
sym.pprint(B)
print ("")

C = 2*b/m

print ("")
print ("C:")
print ("")
sym.pprint(C)
print ("")

A, B, C = sym.symbols('A, B, C', real=True, positive=True)
s, t = sym.symbols('s, t')

transfer_function = A/((s**2) + (C*s) + B)

print ("")
print ("Transfer function:")
print ("")
sym.pprint(transfer_function)
print ("")

# preparing timespan for graphing

n_points = 500
t_final = 0.2
t_span = np.linspace(0, t_final, n_points)

# impulse response

dx2_imp_x = 1
imp_x = transfer_function * dx2_imp_x

print ("")
print ("imp x")
sym.pprint(imp_x)
print ("")

t_imp_x = sym.inverse_laplace_transform(imp_x, s, t)

print ("")
print ("Impulse response for Transfer function")
print ("")
sym.pprint(t_imp_x)
print ("")

# step response

step_x = 1 / s
step_response_x = transfer_function  * step_x

print ("")
print ("step x")
sym.pprint(step_response_x)
print ("")

t_step_x = sym.inverse_laplace_transform(step_response_x, s, t)

print ("")
print ("Step response for Transfer function")
print ("")
sym.pprint(t_step_x)
print ("")



#NOTE:  Need to intialise transfer functions numerically before can be outputted graphically
#CURRENT ISSUE transfer function inputs correctly, however the inverse
#laplace transform results in an infinite recurrance, step and impulse 
#response cant be calcualted
