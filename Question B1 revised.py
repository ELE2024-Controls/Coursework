import sympy as sym
import numpy as np
from scipy.integrate import solve_ivp

m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, x, xe, x1, x2, k1, k2, ieq, y, veq = sym.symbols('m, g, d, δ, r, R, l, L0, L1, α, c, k, b, φ, v, x, xe, x1, x2, k1, k2, ieq, y, veq')

#Question B1 - Show that the system can equilibrate only at those positions xe that satisfy xmin < xe < xmax,

# GIVEN VALUES!

#Ball values
m_value = 0.425
r_value = 0.125
g_value = 9.81

#Measurements
d_value = 0.42
δ_value = 0.65

# Magnet values
R_value = 53
L0_value = 0.120
L1_value = 0.025
α_value = 1.2
c_value = 6815

# Dampner/spring values
k_value = 1880
b_value = 10.4

# Board values
φ_value = 0.733038 #radians

#END OF INITIALISATION


#== SIMULATION FRAMEWORKS =====
#Will create library later on for this for other questions
#Find non-linear and linear equations
#For linear simulation, use transfer functions.
class System: #For Non-linear simulation

    def __init__(self, mass=m_value, radius=r_value,dampening=b_value, spring=k_value, x=δ_value/2, velocity=0): #initialising x between magnet and damp.
        """
        Initialises ball class and it's parameters.
        :param mass: mass of ball (kg)
        :param velocity: Initial velocity (m/s)
        :param x: Initial x-position (m)
        :param radius: size of radius in ball (m)
        :param dampening : dampener coefficient(Ns/m)
        :param spring : spring coefficient(N/m)
        """
        self.mass = mass
        self.radius = radius
        self.dampening = dampening
        self.spring = spring
        self.x = x
        self.velocity = velocity

    def move(self, voltage, dt):
        """
        This function computes and updates new position
        of ball and apply given voltage for "dt" time
        :param voltage: voltage input (V)
        :param dt: Time discrete interval
        """

        def system_dynamics(t, z): #TBD, insert dynamics for IVP solve. Placeholder below
             theta = z[2]
            # return [self.velocity * np.cos(theta),
            #         self.velocity * np.sin(theta),
            #         self.velocity * np.tan(steering_angle_rad) / self.length]

        # next we need to solve the IVP
        z_initial = [self.x, self.y, self.theta]
        solution = solve_ivp(system_dynamics,
                             [0, dt],
                             z_initial)
        self.x = solution.y[0][-1]
        self.y = solution.y[1][-1]
        self.theta = solution.y[2][-1]

# ========== END OF SIMULATION FRAMEWORK ===========

# Equilibrium point
# x - position of ball
# by graph of example, ball starts at, if δ max length
x2eq = 0

# diff x2 = (v**2)*(2*c/m*3*(R*y+L)**2) + (2/3)*g*sym.sin(φ) - (2/3)*b*x2 - (2/3*m)*k1*(x1 - d)-(2/3*m)*k2*(x1-d)**3 = 0 at equilibrium point

y = δ - x1
L = L0 + (L1 * sym.exp(-α * y))

x1 = x
x2 = sym.diff(x1)
dx2 = sym.diff(x2)
dx2 = 2*c*v**2/(R*y+L)**2 + 2*g*sym.sin(φ) - 2*k*x1/m - 2*b*x2/m

# defining Xe(x1) and X2 at Equilibrium point

x1eq = xe
x2eq = 0

dx2_eq = 2*c*veq**2/(R*y+L)**2 + 2*g*sym.sin(φ) - 2*k*x1eq/m - 2*b*x2eq/m

print("")
print("ẋ2 at equilibrium =")
print("")
sym.pprint(dx2_eq)

dx2_eq = 0
X_min = d + (m * g * sym.sin(φ)/k)
X_max = δ

# xmin < xe < δ,

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

#finding values of L and Y to be used in the determined functions

y_min_value = δ_value - X_min_value
y_max_value = δ_value - X_max_value

L_min_value = L0_value + (L1_value * sym.exp(-α_value * y_min_value))
L_max_value = L0_value + (L1_value * sym.exp(-α_value * y_max_value))

#checking validity of dx2eq at xmin and max

dx2eq_min = 2*c_value*veq**2/(R_value*y_min_value+L_min_value)**2 + 2*g_value*sym.sin(φ_value) - 2*k_value*X_min_value/m_value - 2*b_value*x2eq/m_value

dx2eq_max = 2*c_value*veq**2/(R_value*y_max_value+L_max_value)**2 + 2*g_value*sym.sin(φ_value) - 2*k_value*X_max_value/m_value - 2*b_value*x2eq/m_value

print("")
print("Validating dx2eq at X_min")
print("")
sym.pprint(dx2eq_min)
print("")
print("Validating dx2eq at X_max")
print("")
sym.pprint(dx2eq_max)

#Determine the equilibrium voltage and current as a function of xe
#and determine the position xe* where the corresponding equilibrium
#voltage attains its maximum value.


#-------------FINDING VEQ AND IEQ IN TERMS OF XE -------------------

y = δ - xe
L = L0 + (L1 * sym.exp(-α * y))

dx2_eq = 2*c*veq**2/(R*y+L)**2 + 2*g*sym.sin(φ) - 2*k*xe/m - 2*b*x2eq/m

veq_eqn = sym.solve(dx2_eq, veq)

print("")
print("Finding veq in terms of xe")
print("")
sym.pprint(veq_eqn)
print("")

veq_eqn = ieq / R

ieq_eqn = sym.solve (veq_eqn, ieq)
print("")
print("Finding Ieq in terms of xe")
print("")
sym.pprint(ieq_eqn)
#-------------------------------------------------------------------------


# CURRENT ERROR  solving gives an expression for Veq but applying v=ir  doesnt work
# TO DO : validation of xmin - xmax
#         xe* where veq is at max value

#maybe do timespan against veq against xe?
# or veq is max when dveq/dxe is 0?
