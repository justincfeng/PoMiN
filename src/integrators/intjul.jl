#-----------------------------------------------------------------------
#   JULIA INTEGRATOR FUNCTION
#-----------------------------------------------------------------------
#
#   This file provides an interface to Julia's OrdinaryDiffEq.jl
#   integrators for use with PoMiN. It supports both autonomous and
#   non-autonomous systems, and handles type conversion properly.
#
#-----------------------------------------------------------------------

using OrdinaryDiffEq

#-----------------------------------------------------------------------

"""
    jlintegrator(z0, F::Function, tspan::Tuple{Real,Real};
                    p=nothing, abstol::Real=1e-10, 
                    reltol::Real=1e-10, integrator=Tsit5())

Integrates ODE using Julia's OrdinaryDiffEq.jl integrators with full
solution object return.

# Arguments
- `z0::RealVec`: Initial state vector (phase space coordinates)
- `F::Function`: Function F(du,u,p,t) defining the ODE system, where:
  - `du`: Output array for derivatives (modified in-place)
  - `u`: Current state vector
  - `p`: Parameters (optional)
  - `t`: Current time
- `tspan::Tuple{Real,Real}`: Time span for integration (t_start, t_end)
- `p`: Optional parameters to pass to the ODE function
- `atol::Real`: Absolute tolerance (default: 1e-10)
- `rtol::Real`: Relative tolerance (default: 1e-10)
- `integrator`: OrdinaryDiffEq.jl integrator to use (default: Tsit5())

# Returns
- `ODESolution`: Complete solution object

# Notes
This function provides access to the full OrdinaryDiffEq.jl solution
object, which includes interpolation capabilities and detailed solver
statistics. For post-Minkowskian N-body problems, the state vector z0
typically contains both position and momentum coordinates for all
particles.
"""
function jlintegrator(z0, F::Function, tspan::Tuple{Real,Real},
                      p=nothing, atol::Real=1e-10, rtol::Real=1e-10,
                      integrator=Tsit5())
    
    # Create and solve ODE problem
    prob = ODEProblem(F, z0, tspan, p)
    return solve(prob, integrator, abstol=atol, reltol=rtol)
end #-------------------------------------------------------------------