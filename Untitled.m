
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

xe = [0.421483921220710:0.01:0.65];
veq = sqrt((-g*sin(phi)+(k*xe/m)).*exp(-2*alpha*delta)/c) .* ((L0*exp(alpha*delta)) + (L1*exp(xe*alpha)) - (R*xe*exp(alpha*delta)) + (R*delta* exp(alpha*delta)) )
plot(xe, veq)

maxveq = max(veq)
