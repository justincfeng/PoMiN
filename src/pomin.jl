#-----------------------------------------------------------------------
#
#   PoMiN: A Hamiltonian Post-Minkowski N-Body Code
#
#       Version 2.0
#
#-----------------------------------------------------------------------

module pomin

using LinearAlgebra
using CSV
using ForwardDiff
using OrdinaryDiffEq

include("pomin-types.jl")   # data type definitions
include("idxer.jl")         # routines for index management
include("HamPM.jl")         # post-Minkowski Hamiltonian
#include("exf.jl")           # External forces
include("integrators/epsi.jl")          # Symplectic Integrator (Tao 2016)
include("integrators/rk4i.jl")          # 4th order Runge-Kutta
include("integrators/int.jl")           # Julia integrator functions
include("gwsc.jl")          # GW strain calculator

end

# test
#pomin.H(3,rand(Float64,2),rand(Float64,12))

masses = [1 0.00000300336937]
Z_init = [0 0 0 0 0 0 101306075.1 0 0 0 0.0000000002983946294 0]
tadap = t -> 1.0  # effectively turns off adaptive time stepping since δ will be multiplied by 1 each timestep
tspan = (0, 64071141430000)
maxit = 1000000
δ = 6407114143
print(pomin.Jsympl(pomin.dH(3,masses,Z_init)))
#sol = pomin.hrkintegrator(Z_init, x -> pomin.dH(3,masses,x) , δ, tadap, tspan, maxit)