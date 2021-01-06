import sympy as sym

m, g, d, δ, r, R, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i, y, Xe = sym.symbols('m, g, d, δ, r, R, L0, L1, α, c, k, b, φ, v, x1, x2, k1, k2, i, y, Xe')

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

# Equilibrium point

Feq = 0
x2eq = 0

X_min = d + (m * g * sym.sin(φ)/k)
X_max = δ

print("")
print("Equation for X_min")
print("")
sym.pprint(X_min)
print("")
print("Equation for X_max")
print("")
sym.pprint(X_max)

X_min_value = d_value + (m_value * g_value * sym.sin(φ_value)/k_value)
X_max_value = δ_value

print("")
print("Calculated X_min")
print("")
sym.pprint(X_min_value)
print("")
print("Calculated X_max")
print("")
sym.pprint(X_max_value)

y_min_value = δ_value - X_min_value
y_max_value = δ_value - X_max_value

L_min_value = L0_value + (L1_value * sym.exp(-α_value * y_min_value))
L_max_value = L0_value + (L1_value * sym.exp(-α_value * y_max_value))

#checking validity of Feq at xmin and max

Feq_min = (v**2)*(2*c_value/m_value*3*(R_value*y_min_value+L_min_value)**2) + (2/3)*g_value*sym.sin(φ_value) - (2/3)*b_value*x2eq - (2/3*m_value)*k1*(X_min_value - d_value)-(2/3*m_value)*k2*(X_min_value-d_value)**3
Feq_max = (v**2)*(2*c_value/m_value*3*(R_value*y_max_value+L_max_value)**2) + (2/3)*g_value*sym.sin(φ_value) - (2/3)*b_value*x2eq - (2/3*m_value)*k1*(X_max_value - d_value)-(2/3*m_value)*k2*(X_max_value-d_value)**3

print("")
print("Validating Feq at X_min")
print("")
sym.pprint(Feq_min)
print("")
print("Validating Feq at X_max")
print("")
sym.pprint(Feq_max)


