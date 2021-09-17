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
include("exf.jl")           # External forces
include("integrators/epsi.jl")          # Symplectic Integrator (Tao 2016)
include("integrators/rk4i.jl")          # 4th order Runge-Kutta
include("integrators/int.jl")           # Julia integrator functions
include("gwsc.jl")          # GW strain calculator

ΦC(rand(Float64,4),rand(Float64),rand(Float64))

end # module
