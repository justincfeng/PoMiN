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
using DoubleFloats

include("pomin-types.jl")   # data type definitions
include("pomin-io.jl")      # input/output routines
include("idxer.jl")         # routines for index management
include("HamPM.jl")         # post-Minkowski Hamiltonian
#include("exf.jl")           # External forces
include("integrators/epsi.jl")          # Symplectic Integrator (Tao 2016)
include("integrators/rk4i.jl")          # 4th order Runge-Kutta
include("integrators/int.jl")           # Julia integrator functions
include("gwsc.jl")          # GW strain calculator

# include("starshot.jl")

# runExperiment(3)

# v_init1 = Double64[-7.00365851051320E-02, -5.35052039083642E-02, -1.70575701379575E-01]  # closest approach ~= 6 AU
# v_init2 = Double64[-0.07295215750065391, -0.05573727243527886, -0.1776832551972042] # closest approach = 4.15 AU
# v_init3 = Double64[-0.07130599062805795, -0.05447701998124462, -0.17367024572512088] # closest approach = 4.71 AU
# v_init4 = Double64[-0.06926018403390836, -0.05291081619979526, -0.16868299845569543] # closest approach = 7.24 AU
# println(v_init2)
# println(secantMethod(closestApproach, v_init2, v_init3, 1, 20))
# println(closestApproach(v_init) * 1.47669196951425 * 6.6845871226706E-09)
# println(stderr, closestApproach(v_init2) * 1.47669196951425 * 6.6845871226706E-09)
# minvec = closestApproach(v_init4) * 1.47669196951425 * 6.6845871226706E-09
# println(norm(minvec))

#println(ForwardDiff.jacobian(x -> closestApproach(x),v_init4))

#println(BroydenMethod(v_init4))

# include("tests/broyden_test.jl")

include("tests/Earth_Moon_SunPotential.jl")

end

