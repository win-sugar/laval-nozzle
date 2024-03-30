
using DataFrames
using CSV
using Plots
using StatsPlots


gasconst = 8.31446262 # m2 kg s-2 K-1 mol-1
MW = 29e-3  # [kg/mol]
cp = 1005  # [J/kg/K]
R = gasconst / MW
gamma = cp/(cp-R)

calculate_mach(gamma,R,T,u) = u/sqrt(gamma*R*T)


# x T p rho ux uy uz
dfc = CSV.read(IOBuffer(replace(read("cal.dat"), UInt8('\t') => UInt8(' '))),
              header=1, delim=" ", ignorerepeated=true, types=Float64, DataFrame)

dfc[!,:M] = calculate_mach.(gamma, R, dfc[!,:T], dfc[!,:ux])

#  area  M  p  T  rho  u
dfa = CSV.read(IOBuffer(replace(read("ana.dat"), UInt8('\t') => UInt8(' '))),
              header=1, delim=" ", ignorerepeated=true, types=Float64, DataFrame)

# temperature

cal = @df dfc plot(:x, [:T],legend=true,xlabel="x[m]",ylabel="T[K]",title="Temperature(pe-0.89)",
         xlims=(0,0.254),lc=:red,lw=2,ls=:dash, size=(600,400), label="cal",
         mode="markers", mc=:red, ms=2, ma=0.5, st="scatter")

@df dfa plot!(cal, :x, [:T],legend=true,xlabel="x[m]",ylabel="T[K]",title="Temperature(pe-0.89)",
         xlims=(0,0.254),lc=:blue,lw=1,ls=:solid, size=(600,400), label="ana" )

savefig("xy-temp.png")

# pressure

cal = @df dfc plot(:x, [:p],legend=true,xlabel="x[m]",ylabel="p[Pa]",title="Pressure(pe-0.89)",
         xlims=(0,0.254),lc=:red,lw=2,ls=:dash, size=(600,400), label="cal",
         mode="markers", mc=:red, ms=2, ma=0.5, st="scatter")

@df dfa plot!(cal, :x, [:p],legend=true,xlabel="x[m]",ylabel="p[Pa]",title="Pressure(pe-0.89)",
         xlims=(0,0.254),lc=:blue,lw=1,ls=:solid, size=(600,400), label="ana" )

savefig("xy-pres.png")

# density

cal = @df dfc plot(:x, [:rho],legend=true,xlabel="x[m]",ylabel="rho[kg/m3]",title="Density(pe-0.89)",
         xlims=(0,0.254),lc=:red,lw=2,ls=:dash, size=(600,400), label="cal",
         mode="markers", mc=:red, ms=2, ma=0.5, st="scatter")

@df dfa plot!(cal, :x, [:rho],legend=true,xlabel="x[m]",ylabel="rho[kg/m3]",title="Density(pe-0.89)",
         xlims=(0,0.254),lc=:blue,lw=1,ls=:solid, size=(600,400), label="ana" )

savefig("xy-dens.png")

# ux

cal = @df dfc plot(:x, [:ux],legend=true,xlabel="x[m]",ylabel="u[m/s]",title="Velocity(pe-0.89)",
         xlims=(0,0.254),lc=:red,lw=2,ls=:dash, size=(600,400), label="cal",
         mode="markers", mc=:red, ms=2, ma=0.5, st="scatter")

@df dfa plot!(cal, :x, [:u],legend=true,xlabel="x[m]",ylabel="u[m/s]",title="Velocity(pe-0.89)",
         xlims=(0,0.254),lc=:blue,lw=1,ls=:solid, size=(600,400), label="ana" )

savefig("xy-velo.png")

# mach number

cal = @df dfc plot(:x, [:M],legend=true,xlabel="x[m]",ylabel="M[-]",title="Mach Number(pe-0.89)",
         xlims=(0,0.254),lc=:red,lw=2,ls=:dash, size=(600,400), label="cal",
         mode="markers", mc=:red, ms=2, ma=0.5, st="scatter")

@df dfa plot!(cal, :x, [:M],legend=true,xlabel="x[m]",ylabel="M[-]",title="Mach Number(pe-0.89)",
         xlims=(0,0.254),lc=:blue,lw=1,ls=:solid, size=(600,400), label="ana" )

savefig("xy-mach.png")

