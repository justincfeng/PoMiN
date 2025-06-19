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

# Core functionality
include("core/pomin-types.jl")                 # Type definitions

# Physics modules
include("physics/Hamiltonians/HamPM.jl")       # Post-Minkowski Hamiltonian
include("physics/gravitational_waves/gwsc.jl")  # GW strain calculator
include("physics/post_minkowski/HO.jl")        # Harmonic oscillator

# Numerical methods
include("integrators/epsi.jl")          # Symplectic Integrator (Tao 2016)
include("integrators/rk4i.jl")          # 4th order Runge-Kutta
include("integrators/intjul.jl")        # Julia integrator functions
include("integrators/tadap.jl")         # Time adaptation functions

# Utilities
include("utils/idxer.jl")               # Index management
include("utils/idgen.jl")               # ID generation
include("utils/broyden.jl")             # Broyden method solver
include("utils/misc.jl")                # Miscellaneous utilities

# Input/Output
include("io/io.jl")                     # I/O routines

# Export public interface

function solve(system::Particles, params::Parameters)
    m = system.m
    z = vcat(system.q...,system.p...)
    
    if params.jli[1] && params.rkl[1] == false 
        include("integrators/intjul.jl")
        return jlintegratorfull(z,(Z,p,t)->HamPM.FHE(Z,p),params.tspan,params=m,abstol=params.jli[2],reltol=params.jli[2],integrator=params.jli[3])
    elseif params.rkl[1] && params.jli[1] == false
        include("integrators/rk4i.jl")
        include("integrators/tadap.jl")
        return rkl_solv(z,HamPM.FHE,m,params.tspan)
    elseif params.jli[1] == false && params.rkl[1] == false
        error("No integrator selected")
    elseif params.jli[1] == true && params.rkl[1] == true
        error("Both integrators selected")
    end
end #-------------------------------------------------------------------

# Types
export RealVec, Particles, Parameters, soln

# Core functionality
export solve

# Integration methods
export jlintegrator, jlintegratorfull

# Initial data generation
export setup_binary_system, setup_circular_orbit
export setup_massless_test_particle, merge_particle_systems
export add_particle, to_com_frame, setup_scattering
export to_phase_tuple

end
