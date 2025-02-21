#-----------------------------------------------------------------------
#   JULIA INTEGRATOR FUNCTIONS
#-----------------------------------------------------------------------
#
#   This file provides a clean interface to Julia's OrdinaryDiffEq.jl
#   integrators for use with PoMiN. It supports both autonomous and
#   non-autonomous systems, and handles type conversion properly.
#
#   Main functions:
#   - jlintegrator: Main integration function, returns solution arrays
#   - jlintegratorfull: Extended version that also returns solver object
#
#   Example usage:
#   ```julia
#   # For autonomous system F(z)
#   z0 = [1.0, 0.0]  # initial conditions
#   function F(du,u,p,t)  # Note: du is modified in-place
#       du[1] = u[2]
#       du[2] = -u[1]  # simple harmonic oscillator
#   end
#   tspan = (0.0, 10.0)
#   sol = jlintegrator(z0, F, tspan)
#
#   # With parameters
#   function F(du,u,p,t)
#       ω = p[1]  # angular frequency
#       du[1] = u[2]
#       du[2] = -ω^2 * u[1]  # harmonic oscillator with frequency ω
#   end
#   sol = jlintegrator(z0, F, tspan, [2.0])  # ω = 2.0
#   ```
#-----------------------------------------------------------------------

using OrdinaryDiffEq

#-----------------------------------------------------------------------

"""
    jlintegratorfull( z0::RealVec , F::Function 
                    , tspan::Tuple{Real,Real} 
                    , params=nothing , abstol::Real=1e-10 
                    , reltol::Real=1e-10 , integrator=Tsit5() )

Integrates ODE using Julia's OrdinaryDiffEq.jl integrators.

Parameters:
- z0: Initial state vector
- F: Function defining the ODE system. Must be of the form F(du,u,p,t) where:
     * du: Output array for derivatives (modified in-place)
     * u: Current state vector
     * p: Parameters (optional)
     * t: Current time
- tspan: Time span for integration (t_start, t_end)
- params: Optional parameters to pass to the ODE function
- abstol: Absolute tolerance for adaptive stepping
- reltol: Relative tolerance for adaptive stepping
- integrator: OrdinaryDiffEq.jl integrator to use (default: Tsit5())

Returns:
- ODESolution: Solution object containing trajectory and metadata
"""
function jlintegratorfull(z0::RealVec, F::Function, tspan::Tuple{Real,Real}
                       , params=nothing, abstol::Real=1e-10
                       , reltol::Real=1e-10, integrator=Tsit5())
    
    # Create and solve ODE problem
    prob = ODEProblem(F, z0, tspan, params)
    return solve(prob, integrator, abstol=abstol, reltol=reltol)
end

#-----------------------------------------------------------------------

"""
    jlintegrator( z0::RealVec , F::Function 
                , tspan::Tuple{Real,Real}
                , params=nothing , abstol::Real=1e-10 
                , reltol::Real=1e-10 , integrator=Tsit5() )

Simplified version of jlintegratorfull that returns the solution.
See jlintegratorfull for full documentation.
"""
function jlintegrator(z0::RealVec, F::Function, tspan::Tuple{Real,Real}
                   , params=nothing, abstol::Real=1e-10
                   , reltol::Real=1e-10, integrator=Tsit5())
    
    return jlintegratorfull(z0, F, tspan, params, abstol, reltol, integrator)
end
