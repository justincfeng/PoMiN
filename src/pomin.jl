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
#using Logging
using DoubleFloats

# Core functionality
include("core/pomin-types.jl")                 # Type definitions

# Physics modules
include("physics/Hamiltonians/HamPM.jl")       # Post-Minkowski Hamiltonian
include("physics/gravitational_waves/gwsc.jl") # GW strain calculator

# Physics utilities
include("physics/Hamiltonians/idxer.jl")       # Phase space indexing
include("physics/initial_data/idgen.jl")       # Initial data generation

# Numerical methods
include("integrators/tadap.jl")                 # Time adaptation functions
include("integrators/rk4i.jl")                  # 4th order Runge-Kutta
include("integrators/intjul.jl")                # Julia integrator functions

# Input/Output
include("io/io.jl")                             # I/O routines

"""
    solve(system::Particles, params::Parameters)

Solve the post-Minkowskian N-body problem for the given particle system and parameters.

# Arguments
- `system::Particles`: The particle system containing masses, positions, and momenta
- `params::Parameters`: Integration parameters and settings

# Returns
- For RK4 integrator: `soln` structure containing the time evolution
- For Julia integrators: `ODESolution` object from OrdinaryDiffEq.jl

# Notes
This is the main entry point for PoMiN simulations. The function automatically
selects the appropriate integrator based on `params.rkl` and constructs the
phase space vector from the particle system.
"""
function solve(system::Particles, params::Parameters)
    m = system.m
    z = vcat(system.q..., system.p...)

    if params.rkl
        # Use RK4 integrator
        return hrkintegrator(params.d, length(m), z, 
                           (Z) -> HamPM.dH(Z, m, params.d),
                           params.δ,
                           (dt, Z, Zdot) -> tadap.tcour(dt, Z, Zdot, params.courant, params.d),
                           params.tspan, params.iter, params.Nrec)
    else
        # Use Julia OrdinaryDiffEq.jl integrator
        return jlintegrator(z, 
                          (du, u, p, t) -> begin
                              du .= HamPM.FHE(u, p)
                          end,
                          params.tspan, m, params.atol, params.rtol,
                          eval(Symbol(params.integrator))(), params.Nrec)
    end
end #-------------------------------------------------------------------

# Types
export RealVec, Particles, Parameters, soln

# Core functionality
export solve

# Integration methods
export jlintegrator, tcour

# Initial data generation
export setup_binary_system, setup_circular_orbit
export setup_massless_test_particle, merge_particle_systems
export add_particle, to_com_frame, setup_scattering

end
