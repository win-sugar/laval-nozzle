
using NLsolve

gasconst = 8.31446262 # m2 kg s-2 K-1 mol-1
MW = 29e-3
cp = 1005
R = gasconst / MW
gamma = cp/(cp-R)

A0 = 0.0016129 # =2.5[inch^2]
p0 = 6895
T0 = 100
rho0 = p0/(R*T0)
a0 = sqrt( gamma*R*T0 )

Ts = 2/(gamma+1) * T0
ps = (2/(gamma+1))^(gamma/(gamma-1)) * p0
rhos = (2/(gamma+1))^(1/(gamma-1)) * rho0
as = sqrt( 2/(gamma+1) ) * a0
us = as
As = 0.00064516 # =1.0[inch^2]
mdot = rhos*us*As

Ae = 0.00096774 # =1.5[inch^2]
pe = 6137

f(x,gamma,p0,pe) = ( p0/pe
                   - ( (1+(gamma-1)/2*x[1]^2) )^(gamma/(gamma-1)) )

g(x) = f(x,gamma,p0,pe)

sol = nlsolve( g, [ 0.3 ])
result = sol.zero
Me = result[1]

Te = T0 / ( (1+(gamma-1)/2*Me^2) )
rhoe = rho0 / ( (1+(gamma-1)/2*Me^2)^(1/(gamma-1)) )
ae = sqrt( R*gamma*Te )
ue = Me*ae

n = 100
for i in  0:n
        x = i*(0.254/n)
        xin = i*(10/n)
        area = ( x<=0.127 ? (1.75-0.75*cos((0.2*xin-1.0)*pi) )*2.54e-2^2
                          : (1.25-0.25*cos((0.2*xin-1.0)*pi) )*2.54e-2^2 )

f(x,area,gamma,Me,Ae) = ( area/Ae
       - ( (1+(gamma-1)/2*x[1]^2)/(1+(gamma-1)/2*Me^2) )^((gamma+1)/(2*(gamma-1)))*Me/x[1] )

g(x) = f(x,area,gamma,Me,Ae)

local sol = nlsolve( g, [ 0.3 ], show_trace=false )

local result = sol.zero
M = result[1]

rho = ( (1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*Me^2) )^(-1/(gamma-1)) * rhoe
T = ( (1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*Me^2) )^(-1) * Te
p = ( (1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*Me^2) )^(-gamma/(gamma-1)) * pe
u = M/Me*( (1+(gamma-1)/2*Me^2)/(1+(gamma-1)/2*M^2) )^(0.5) * ue

println("$x  $area  $M  $p  $T  $rho  $u")

end



