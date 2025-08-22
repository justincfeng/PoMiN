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

∂ = (f,Z)->ForwardDiff.gradient(f,Z)

# Core functionality
include("core/pomin-types.jl")                 # Type definitions

# Physics modules
include("core/physics/Hamiltonians/HamPM.jl")    # Post-Minkowskian Hamiltonian
include("core/physics/Hamiltonians/HamTools.jl") # Hamiltonian tools
include("post/gravitational_waves/gwsc.jl") # GW strain calculator

# Physics utilities
include("core/initial_data/idgen.jl")    # Initial data generation

# Integrators
include("core/integrators/rk4i.jl")              # RK4 integrator
include("core/integrators/intjul.jl")            # Julia ODE integrators

# Input/Output
include("utils/io.jl")                             # I/O routines

"""
    solve(system::Particles, params::Parameters,
          Phi::Function=z->zero(typeof(z[1])))

Solve the post-Minkowskian N-body problem for the given particle system and parameters.

# Arguments
- `system::Particles`: Specifies masses, positions, and momenta
- `params::Parameters`: Integration parameters and settings
- `Phi::Function`: Specifies an external potential function

# Returns
- For RK4 integrator: `soln` structure containing the time evolution
- For Julia integrators: `ODESolution` object from OrdinaryDiffEq.jl

# Notes
This is the main entry point for PoMiN simulations. The function
automatically selects the appropriate integrator based on `params.rkl`
and constructs the phase space vector from the particle system.
"""
function solve(system::Particles, params::Parameters, 
               Phi::Function=z->zero(typeof(z[1])))
    m = system.m
    z = vcat(system.q..., system.p...)

    if params.rkl
        # Use RK4 integrator
        return hrkintegrator(params.d, length(m), z, 
                           (Z) -> HamPM.dH(Z, m, params.d)+∂(Phi,Z),
                           params.δ,
                           params.tspan, params.iter,
                           (dt, Z, Zdot) -> tcour(dt, Z, Zdot, params.courant, params.d),
                           params.Nrec)
    else
        # Use Julia OrdinaryDiffEq.jl integrator
        return jlintegrator(z, 
                          (du, u, p, t) -> begin
                              du .= HamPM.FHE(u, m, params.d) + 
                                    HamPM.Jsympl(∂(Phi,u))
                          end,
                          params.tspan, m, params.atol, params.rtol,
                          eval(Symbol(params.integrator))(), params.Nrec)
    end
end #-------------------------------------------------------------------

# Types
export RealVec, Particles, Parameters, soln

# Core functionality
export solve

# Physics modules
export HamPM, dp_scatter

# Integration methods
export jlintegrator, hrkintegrator, tcour, rk4map, tnone

# Parameter constructors
export Parameters, ParametersRK4, ParametersJulia

# Initial data generation
export setup_binary_system, setup_circular_orbit, setup_elliptical_orbit
export setup_massless_test_particle, merge_particle_systems
export add_particle, to_com_frame, setup_scattering

end
