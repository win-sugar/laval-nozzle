
using NLsolve

gasconst = 8.31446262 # m2 kg s-2 K-1 mol-1
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
mdot = rhos*us*As

Ae = 0.00096774
pe = 5171
Me = 0.50200727578001
ae = gamma*pe*Me*Ae/mdot
Te = ae^2/(gamma*R)
ue = Me * ae
rhoe = pe/(R*Te)

n = 100
for i in  75:n
    if i == 75  # <-- shock position
        area = 0.0008129027418810686
        x = 0.1921177321640935
    else
        x = i*(0.254/n)
        xin = i*(10/n)
        area = ( x<=0.127 ? (1.75-0.75*cos((0.2*xin-1)*pi) )*2.54e-2^2
                          : (1.25-0.25*cos((0.2*xin-1)*pi) )*2.54e-2^2 )
    end

f(x,area,gamma,Me,Ae) = ( area/Ae
      - ( (1+(gamma-1)/2*x[1]^2)/(1+(gamma-1)/2*Me^2) )^((gamma+1)/(2*(gamma-1)))*Me/x[1] )

g(x) = f(x,area,gamma,Me,Ae)

sol = nlsolve( g, [ 0.3 ], show_trace=false, ftol=1e-15)

result = sol.zero
M = result[1]

rho = ( (1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*Me^2) )^(-1/(gamma-1)) * rhoe
T = ( (1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*Me^2) )^(-1) * Te
p = ( (1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*Me^2) )^(-gamma/(gamma-1)) * pe
u = M/Me*( (1+(gamma-1)/2*Me^2)/(1+(gamma-1)/2*M^2) )^(0.5) * ue

println("$x  $area  $M  $p  $T  $rho  $u")

end



