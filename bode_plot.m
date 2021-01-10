
m = 0.425
r=0.125
g=9.81
d=0.42
delta=0.65
R=53
L0=0.12
L1=0.025
alpha = 1.2
c=6815
k = 1880
b=10.4
phi = 0.733838


xmin = d+(m*g*sin(phi)/k)

xmax = delta

xe = 0.75*xmin + 0.25*xmax

y= delta - xe
 
L=L0+(L1*exp(-delta*y))

veq = sqrt((-g*sin(phi)+(k*xe/m))*exp(-2*alpha*delta)/c) * ((L0*exp(alpha*delta)) + (L1*exp(xe*alpha)) - (R*xe*exp(alpha*delta)) + (R*delta* exp(alpha*delta)) )

A= 2*c/m*(R*y+L)^2*(2*veq)

B= 2*k/m

C=2*b/m

sys = tf([A],[1, C, B])


bode(sys)
grid