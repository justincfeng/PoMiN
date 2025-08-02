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

"""    jlintegrator(z0, F::Function, tspan::Tuple{Real,Real};
                    p=nothing, abstol::Real=1e-10, 
                    reltol::Real=1e-10, integrator=Tsit5(), Nrec::Int=100)

Integrates ODE using Julia's OrdinaryDiffEq.jl integrators with configurable
solution saving.

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
- `Nrec::Int`: Recording frequency (Nrec=-1 saves only last point, Nrec=-N saves last N points)

# Returns
- `ODESolution`: Solution object (full or truncated based on Nrec)

# Notes
This function provides access to OrdinaryDiffEq.jl solution objects with
configurable saving. When Nrec=-1, only the final state is saved for memory
efficiency in scattering calculations. For post-Minkowskian N-body problems,
the state vector z0 typically contains both position and momentum coordinates.
"""
function jlintegrator(z0, F::Function, tspan::Tuple{Real,Real},
                      p=nothing, atol::Real=1e-10, rtol::Real=1e-10,
                      integrator=Tsit5(), Nrec::Int=100)
    
    # Create and solve ODE problem
    prob = ODEProblem(F, z0, tspan, p)
    
    if Nrec == -1
        # Save only the last point for memory efficiency
        sol = OrdinaryDiffEq.solve(prob, integrator, abstol=atol, reltol=rtol, save_everystep=false)
    elseif Nrec < 0
        # Save last |Nrec| points - solve with full solution then truncate
        sol_full = OrdinaryDiffEq.solve(prob, integrator, abstol=atol, reltol=rtol)
        n_keep = min(abs(Nrec), length(sol_full.t))
        # Create truncated solution with last n_keep points
        t_trunc = sol_full.t[end-n_keep+1:end]
        u_trunc = sol_full.u[end-n_keep+1:end]
        sol = remake(sol_full, t=t_trunc, u=u_trunc)
    else
        # Default behavior: save full solution or with specified sampling
        sol = OrdinaryDiffEq.solve(prob, integrator, abstol=atol, reltol=rtol)
    end
    
    return sol
end #-------------------------------------------------------------------