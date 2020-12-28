import sympy as sym
import numpy as np

m, g, d, δ, r, R, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i = sym.symbols('m, g, d, δ, r, R, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i')

L = L0 + (L1 * exp(-α * y))
Fmag = c(i**2/y**2)
y = δ - x1

k = k1
k2 = F.diff(k1)
x2 = F.diff(x1)

F = ((veq^2)*(2*c/(3*(R*y+L)**2*m))) + ((2/(3*g))*np.sin(φ)) - ((2/3)*b*x2eq) - ((2/(3*m)*(k1*(x1eq - d)))-((2/(3*m)*k2*(x1eq-d)**3))

X_min = d + ((m * g * (sym.sin(φ))/k)
X_max = δ


# Equilibrium point

Feq = 0
x2eq = 0
x1eq = Xe
veq = ieq * R


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
    return float(z.subs([(m, m_value), (g, g_value), (d, d_value), (δ, δ_value), (r, r_value), (R, R_value), (L0, L0_value), (L1, L1_value), (α, α_value), (c, c_value), (k, k_value), (b, b_value), (φ, φ_value)]))


X_min_value = evaluate_at_given_parameters(X_min)
X_max_value = evaluate_at_given_parameters(X_max)
F_value_1 = evaluate_at_given_parameters(F, Xe = X_min_value)
F_value_2 = evaluate_at_given_parameters(F, Xe = X_max_value)

sym.pprint(F_value_1)
print("")
sym.pprint(F_value_2)
print("")
sym.pprint(X_min_value)
print("")
sym.pprint(X_max_value)
print("")
