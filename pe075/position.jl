
using NLsolve

gasconst = 8.31446262
MW = 29e-3
cp = 1005
R = gasconst / MW
gamma = cp/(cp-R)

A0 = 0.0016129
p0 = 6895
T0 = 100
rho0 = p0/(R*T0)
a0 = sqrt( gamma*R*T0 )

Ts = 2/(gamma+1) * T0
ps = (2/(gamma+1))^(gamma/(gamma-1)) * p0
rhos = (2/(gamma+1))^(1/(gamma-1)) * rho0
as = sqrt( 2/(gamma+1) ) * a0
us = as
As = 0.00064516

area = 0.0008129027418810686

f(x,area) = area - ( 1.25 - 0.25*cos( (0.2*(x[1]/2.54e-2)-1.0)*pi ) )*2.54e-2^2

g(x) = f(x,area)

ini = 0.18
sol = nlsolve( g, [ ini ], show_trace=true, ftol=1e-15, method=:newton)

result = sol.zero
x = result[1]
println("position = $x")


