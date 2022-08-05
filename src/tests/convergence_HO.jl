using LinearAlgebra
using DoubleFloats

include("HO.jl")

const tpfl = Double64

z0 = tpfl[1; 0]

tspan = (tpfl(0), tpfl(2 * pi / 100))

Δt = tspan[2] - tspan[1]

ω = Double64(64)

nn = 20

q = zeros(tpfl, 9)

Q = zeros(tpfl, 7)

for i = 0:8
    n = nn * 2^i
    δ = tpfl(Δt / n)
    sol = hsintegrator(3, 1, z0, HamHO.dH, δ, ω, tspan, n)
    q[i+1] = sol.z[n][1]
end

for i = 1:7
    Q[i] = abs((q[i] - q[i+1]) / (q[i+1] - q[i+2]))
end

println(Q)
