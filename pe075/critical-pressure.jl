
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

As = 0.00064516
Ts = 2/(gamma+1) * T0
ps = (2/(gamma+1))^(gamma/(gamma-1)) * p0
rhos = (2/(gamma+1))^(1/(gamma-1)) * rho0
as = sqrt( 2/(gamma+1) ) * a0

Ae = 0.00096774

#####################################################

using NLsolve

pe(M) = ps * ( ((gamma+1)/2)/(1+(gamma-1)/2*M^2) )^(gamma/(gamma-1))

function f(x,gamma,As,Ae)
  [ Ae/As - ( (1+(gamma-1)/2*x[1]^2)/((gamma+1)/2) )^((gamma+1)/(2*(gamma-1)))/x[1] ]
end

g(x) = f(x,gamma,As,Ae)

# --> subsonic flow

sol = nlsolve(g, [  0.5 ], show_trace=true)
Me = sol.zero[1]

println("subsonic-Me = ",Me)
println("subsonic-pe = ",pe(Me))

# --> supersonic flow

sol = nlsolve(g, [  1.5 ], show_trace=true)
Me = sol.zero[1]

println("supersonic-Me = ",Me)
println("supersonic-pe = ",pe(Me))


