import sympy as sym
import numpy as np
import matplotlib.pyplot as plt


m, g, d, δ, r, R, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i, y = sym.symbols('m, g, d, δ, r, R, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i, y ')

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

X_min = d_value + (m_value * g_value * sym.sin(φ_value)/k_value)
X_max = δ_value
x_value = 0.75 * X_min + 0.25 * X_max

y_value = δ_value - x_value
L_value = L0_value + (L1_value * sym.exp(-α_value * y_value))

q, u, w = sym.symbols('q, u, w', real=True, positive=True)
s, t = sym.symbols('s, t')

q = 2*c_value/3*m_value(R_value*y_value+L_value)**2

u = 2*b_value/3*m_value

w = 2(k1+k2*3(x_value-d_value)**2)/3*m_value

transfer_function = q / (s + (u/2) - ((u**2 - u*w)**0.5/2))(s + (u/2) - ((u**2 - u*w)**0.5/2))

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

