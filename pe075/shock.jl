
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
pe = 5171

#####################################################

using NLsolve

function f(x,gamma,As,Ae,p0,pe,ps)

# 1  2  3  4  5  6
# M1 M2 A1 p1 p2 Me

  [ x[3]/As - ( (1+(gamma-1)/2*x[1]^2)/((gamma+1)/2) )^((gamma+1)/(2*(gamma-1)))/x[1],
    x[4]/ps - ( ((gamma+1)/2)/(1+(gamma-1)/2*x[1]^2) )^(gamma/(gamma-1)),
    x[5]/x[4] - 1 - (2*gamma/(gamma+1)*(x[1]^2-1)),
    x[2]^2 - ( (1+(gamma-1)/2*x[1]^2)/( gamma*x[1]^2-(gamma-1)/2 ) ),
    x[5]/pe - ( (1+(gamma-1)/2*x[6]^2)/(1+(gamma-1)/2*x[2]^2) )^(gamma/(gamma-1)),
    x[3]/Ae - ( (1+(gamma-1)/2*x[2]^2)/(1+(gamma-1)/2*x[6]^2) )^((gamma+1)/(2*(gamma-1)))*x[6]/x[2] ]
end

g(x) = f(x,gamma,As,Ae,p0,pe,ps)

#                   M1   M2    A1     p1    p2     Me
sol = nlsolve(g, [  2,  0.5,  7e-4,  3e3,  6e3,  0.3 ], show_trace=true)
result = sol.zero

println("M1 = ",result[1])
println("M2 = ",result[2])
println("A1 = ",result[3])
println("p1 = ",result[4])
println("p2 = ",result[5])
println("Me = ",result[6])


