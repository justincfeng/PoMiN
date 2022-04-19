using DoubleFloats
using Plots
using BenchmarkTools

include("pomin-types.jl")
include("integrators/epsi.jl")
include("integrators/rk4i.jl")

include("HO.jl")

z0 = Double64[ 1 ; 0 ]

δ = Double64(0.0001)

tspan = (Double64(0),Double64(2*pi))

maxit = 1000000

ω = Double64(1/2)

tadap = t -> 1.  # effectively turns off adaptive time stepping since δ will be multiplied by 1 each timestep

sol = hrkintegrator( z0, HamHO.dH , δ , tadap, tspan , maxit )
sols = hsintegrator( z0, HamHO.dH , δ , ω , tspan , maxit )

n = length(sol.t)

HV = zeros(Double64,n)
V = zeros(Double64,n)
dV = zeros(Double64,n)
HVs = zeros(Double64,n)
Vs = zeros(Double64,n)
dVs = zeros(Double64,n)

HV0 = HamHO.H(z0)

for i=1:n
    V[i] = sol.z[i][1]
    dV[i] = V[i] - cos(Float64(sol.t[i]))
    HV[i] = HamHO.H(sol.z[i]) - HV0
    Vs[i] = sols.z[i][1]
    dVs[i] = Vs[i] - cos(Float64(sols.t[i]))
    HVs[i] = HamHO.H(sols.z[i]) - HV0
end

plot(sol.t,HVs)
plot(sol.t,V)

plot(sol.t,dV)
plot(sol.t,dVs)
