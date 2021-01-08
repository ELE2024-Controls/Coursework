import  numpy  as np
import  matplotlib.pyplot  as plt
import  scipy.integrate  as spi

def  horizontal_mass_spring(t, y, mass , stiffness , force):
    """ Dynamics  of  horizontal  mass -spring  system

    :param t: time
    :param y: state  vector  where y[0] is the  displacement  ofthe  spring  and y[1] is its  velocity
    :param  mass: mass of  suspended  object:param  stiffness: spring  stiffness
    :param  force: external  force
    :return: right  hand  side of  system  dynamics"""

    return[y[1], (force  - stiffness*y[0])/ mass]

# System parameters
m = 2# mass (in kg)
k = 1000# stiffness (in N/m)
F = 0# external force
t0 = 0# initial time
tf = 10 # final time
y0 = [0.1, 2]# initial condition

num_time_points = 1000
t_span = np.linspace(t0, tf, num_time_points)

# Call the solver
sol = spi.solve_ivp(lambda t, y: horizontal_mass_spring(t, y, m, k, F),[t0 , tf], y0, t_eval=t_span)

plt.figure(1)
plt.plot(sol.t, sol.y[0])
plt.plot(sol.t, sol.y[1])
plt.xlabel('Time(s)')
plt.ylabel('Solution')
plt.legend(['Position','Velocity'])
plt.show()