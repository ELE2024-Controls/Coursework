import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, xe, x1, x2, k1, k2, i, y, veq = sym.symbols('m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, xe, x1, x2, k1, k2, i, y, veq')

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

y = δ - x1
L = L0 + (L1 * sym.exp(-α * y))

x1 = xe
x2 = sym.diff(x1)
dx2 = sym.diff(x2)
dx2 = (v**2)*(2*c/m*3*(R*y+L)**2) + (2/3)*g*sym.sin(φ) - (2/3)*b*x2 - (2/3*m)*k1*(x1 - d)-(2/3*m)*k2*(x1-d)**3

# defining Xe(x1) and X2 at Equilibrium point

x2eq = 0
dx2eq = 0
X_min = d + (m * g * sym.sin(φ)/k)
X_max = δ


X_min = d_value + (m_value * g_value * sym.sin(φ_value)/k_value)
X_max = δ_value
xe = 0.75 * X_min + 0.25 * X_max

y_value = δ_value - xe
L_value = L0_value + (L1_value * sym.exp(-α_value * y_value))

q, u, w = sym.symbols('q, u, w', real=True, positive=True)
s, t = sym.symbols('s, t')

# finding q, u and w, laying out transfer function

q = 2*c_value/3*m_value(R_value*y_value+L_value)**2

u = 2*b_value/3*m_value

w = 2(k1+k2*3(xe-d_value)**2)/3*m_value

transfer_function = q / (s + (u/2) - ((u**2 - u*w)**0.5/2))(s + (u/2) - ((u**2 - u*w)**0.5/2))

# preparing timespan for graphing

n_points = 500
t_final = 0.2
t_span = np.linspace(0, t_final, n_points)

# impulse response

F_imp_x = 1
imp_x = transfer_function  * F_imp_x
t_imp_x = sym.inverse_laplace_transform(imp_x, s, t)

print ("")
print ("Impulse response for Transfer function")
print ("")
sym.pprint(t_imp_x.simplify())
print ("")

plt.plot (t_span, t_imp_x)
plt.suptitle('Impulse response', fontsize=16)
plt.xlabel('time(s)')
plt.ylabel('Position of wooden ball')
plt.grid()
plt.show()

# step response

F_step_x = 1 / s
step_x = transfer_function  * F_step_x
t_step_x = sym.inverse_laplace_transform(step_x, s, t)

print ("")
print ("Step response for Transfer function")
print ("")
sym.pprint(t_step_x.simplify())
print ("")

plt.plot(t_span, t_imp_x)
plt.suptitle('Step response', fontsize=16)
plt.xlabel('time(s)')
plt.ylabel('Position of wooden ball')
plt.grid()
plt.show()


#CURRENT ERROR

# C:/Users/david/PycharmProjects/control coursework/Question B3.py:37: SyntaxWarning: 'int' object is not callable; perhaps you missed a comma?
#  w = 2(k1+k2*3(x_value-d_value)**2)/3*m_value
# C:/Users/david/PycharmProjects/control coursework/Question B3.py:37: SyntaxWarning: 'int' object is not callable; perhaps you missed a comma?
#  w = 2(k1+k2*3(x_value-d_value)**2)/3*m_value
# Traceback (most recent call last):
#  File "C:/Users/david/PycharmProjects/control coursework/Question B3.py", line 33, in <module>
#    q = 2*c_value/3*m_value(R_value*y_value+L_value)**2
# TypeError: 'float' object is not callable
