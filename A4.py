import sympy as sym

m, g, c, V, R, L, theta, k, b, x, x2 = \
    sym.symbols('m, g, c, V, R, L, theta, k, b, x, x2')
print("===========Input eq")
z = (2*c*(V**2))/(m*(4*R+L)**2) + 2*g * sym.sin(theta) - 2*k/m*(x - b) - 2*b/m * x2
sym.pprint(z)

equilibrium_state = ([(x2, 0),(V, 0)])

print("=====Differentiated z wrt x2======")
diff_z_x2 = z.diff(x2)
sym.pprint(diff_z_x2)
print("=====Differentiated z wrt V======")
diff_z_V = z.diff(V)
sym.pprint(diff_z_V)
print("====Substituting for equilibrium values")
print("=====Substituting for Veq======")
diff_z_V_eq = diff_z_V.subs(equilibrium_state)
sym.pprint(diff_z_V_eq)
print("=====Substituting for x2======")
diff_z_x2_eq = diff_z_x2.subs(equilibrium_state)
sym.pprint(diff_z_x2_eq)
