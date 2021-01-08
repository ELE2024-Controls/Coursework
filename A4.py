import sympy as sym

m, g, c, V, R, L, theta, k, b, x, x2 = \
    sym.symbols('m, g, c, V, R, L, theta, k, b, x, x2')
print("===========Input eq")
z = ((2 * m * c * V**2)/(4 * R + L)**2) + 2*m**2*g*sym.sin(theta) - 2*m*k*x - 2*m*b*x2
sym.pprint(z)

def evaluate_at_equilibrium(f):
    f.subs([(x2, 0),(V, 0)])

print("=====Differentiated z wrt x2======")
diff_z_x2 = z.diff(x2)
sym.pprint(diff_z_x2)
print("===========")
diff_z_V = z.diff(V)
sym.pprint(diff_z_V)
print("===========")
eq_z = evaluate_at_equilibrium(z.diff(V))

sym.pprint(eq_z)