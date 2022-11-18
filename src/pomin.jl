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
using Logging

include("pomin-types.jl")   # data type definitions
include("pomin-io.jl")      # input/output routines
include("idxer.jl")         # routines for index management
include("HamPM.jl")         # post-Minkowski Hamiltonian
#include("exf.jl")           # External forces
include("integrators/epsi.jl")          # Symplectic Integrator (Tao 2016)
include("integrators/rk4i.jl")          # 4th order Runge-Kutta
include("integrators/int.jl")           # Julia integrator functions
include("gwsc.jl")          # GW strain calculator

include("tests/AlphaAB_orbit_test.jl")

end

