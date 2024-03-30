
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

n = 100
for i in 0:n
        x = i*(0.254/n)
        xin = i*(10/n)
        area = ( x<=0.127 ? (1.75-0.75*cos((0.2*xin-1)*pi) )*2.54e-2^2
                          : (1.25-0.25*cos((0.2*xin-1)*pi) )*2.54e-2^2 )

f(x,area,gamma,As) = ( area/As
       - ( (1+(gamma-1)/2*x[1]^2)/((gamma+1)/2) )^((gamma+1)/(2*(gamma-1)))/x[1] )

g(x) = f(x,area,gamma,As)

ini = x<=0.127 ? 0.3 : 1.1
sol = nlsolve( g, [ ini ], show_trace=false, ftol=1e-15)
result = sol.zero
M = result[1]

rho = ( ((gamma+1)/2)/(1+(gamma-1)/2*M^2) )^(1/(gamma-1)) * rhos
T = ((gamma+1)/2)/(1+(gamma-1)/2*M^2) * Ts
p = ( ((gamma+1)/2)/(1+(gamma-1)/2*M^2) )^(gamma/(gamma-1)) * ps
u = M*( ((gamma+1)/2)/(1+(gamma-1)/2*M^2) )^0.5 * us

println("$x  $area  $M  $p  $T  $rho  $u")

end



