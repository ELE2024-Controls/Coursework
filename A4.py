import sympy as sym

m, g, c, v, R, L, φ, k, b, x1, x2, y = \
    sym.symbols('m, g, c, V, R, L, theta, k, b, x, x2, y')
print("===========Input eq")
z = 2*c*v**2/(R*y+L)**2 + 2*g*sym.sin(φ) - 2*k*x1/m - 2*b*x2/m
sym.pprint(z)

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
a = diff_z_v.subs(equilibrium_state)
b = diff_z_x1.subs(equilibrium_state)
c = diff_z_x2.subs(equilibrium_state)
sym.pprint(a)
sym.pprint(b)
sym.pprint(c)

def subsitute(z): #Substitution function
    subsititons = [(M, M_value), (m, m_value), (ell, ell_value), (g, g_value)]
    return float(z.subs(subsititons))