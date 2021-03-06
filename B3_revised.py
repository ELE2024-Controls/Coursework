import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import control as Ctrl

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


dx2 = 2*c*v**2/(R*y+L)**2 + 2*g*sym.sin(φ) - 2*k*x1/m - 2*b*x2/m
z = dx2
# defining Xe(x1) and X2 at Equilibrium point

equilibrium_state = ([(x1, 0),(x2, 0),(v, 0)])

print("=====Differentiated z wrt x1======")
diff_z_x1 = z.diff(x1)
sym.pprint(diff_z_x1)
print("=====Differentiated z wrt x2======")
diff_z_x2 = z.diff(x2)
sym.pprint(diff_z_x2)
print("=====Differentiated z wrt V======")
diff_z_v = z.diff(v)
sym.pprint(diff_z_v)
print("====Substituting for equilibrium values")
print("=====Substituting for Veq======")
diff_z_v_eq = diff_z_v.subs(equilibrium_state)
sym.pprint(diff_z_v_eq)
print("=====Substituting for x2======")
diff_z_x2_eq = diff_z_x2.subs(equilibrium_state)
sym.pprint(diff_z_x2_eq)
print("===========")

A = diff_z_v.subs(equilibrium_state)
B = diff_z_x1.subs(equilibrium_state)
C = diff_z_x2.subs(equilibrium_state)

sym.pprint(A)
sym.pprint(B)
sym.pprint(C)


def substiute(z): #Substitution function
    subsititons = [(k, k_value), (m, m_value), (b,b_value)]
    return float(z.subs(subsititons))

A_value = substiute(A)
B_value = substiute(B)
C_value = substiute(C)

sym.pprint(A_value)
sym.pprint(B_value)
sym.pprint(C_value)

a, b, c = sym.symbols('a, b, c', real=True, positive=True, imaginary=False)
s, t = sym.symbols('s, t')
transfer_function = a/((s**2) + (c*s) + b)

print ("")
print ("Transfer function:")
print ("")
sym.pprint(transfer_function)
print ("")

print ("Impulse response")
x1_t = sym.inverse_laplace_transform(transfer_function, s, t)
sym.pprint(x1_t.simplify())
print(sym.latex(x1_t.simplify()))

num_x = [1]
den_x = [1 ,C_value, B_value] #Figure priming transfer function in program
G_x = Ctrl.TransferFunction(num_x, den_x)

t_span = np.linspace(0, 20, 2500)
F_input = np.sin(100*t_span**2)

#No controller simulation
t_out_x, x1_out, _ = Ctrl.forced_response(G_x, t_span, F_input)

plt.plot(t_out_x,x1_out)
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.show()
# x1_freq_t = sym.inverse_laplace_transform(transfer_function*w**2/(s**2 + w**2), s, t)
# sym.pprint(x1_freq_t)
# print(sym.latex(x1_freq_t.simplify())) Takes too long to calculate





