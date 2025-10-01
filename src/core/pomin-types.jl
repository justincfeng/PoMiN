#-----------------------------------------------------------------------
#   TYPE ALIASES
#-----------------------------------------------------------------------

"""
    RealVec{T<:Real}

Type alias for a vector of real numbers. Equivalent to `Array{T}` where
`T` is a subtype of `Real`. Used throughout PoMiN for representing
position and momentum vectors in phase space.
"""
const RealVec{T<:Real} = Array{T}

#-----------------------------------------------------------------------
#   TYPE-CONVERSION FUNCTION
#-----------------------------------------------------------------------
"""
    tpnum(tpfl::Type) -> Tuple

Generate a tuple of numeric constants (0 through 9) for a given floating point type.

# Arguments
- `tpfl::Type`: The floating point type to generate constants for

# Returns
- `Tuple`: A tuple containing constants 0-9 in the specified type

# Examples
```julia
# Generate Float64 constants
constants = tpnum(Float64)  # Returns (0.0, 1.0, 2.0, ..., 9.0)

# Generate BigFloat constants
constants = tpnum(BigFloat)  # Returns (0, 1, 2, ..., 9) as BigFloat
```
"""
function tpnum(tpfl::Type)
    if tpfl <: Real
        return (tpfl(0),tpfl(1),tpfl(2),tpfl(3),tpfl(4),
                tpfl(5),tpfl(6),tpfl(7),tpfl(8),tpfl(9))
    else
        # Use the same precision type instead of defaulting to Float64
        return (tpfl(0),tpfl(1),tpfl(2),tpfl(3),tpfl(4),
                tpfl(5),tpfl(6),tpfl(7),tpfl(8),tpfl(9))
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   PARTICLE STRUCT
#-----------------------------------------------------------------------
"""
    Particles

Datatype representing a collection of particles in the N-body system.

# Fields
- `m::RealVec`: Vector of particle masses
- `q::Array{RealVec}`: Array of position vectors for each particle
- `p::Array{RealVec}`: Array of momentum vectors for each particle

"""
struct Particles
    m::RealVec
    q::Array{RealVec}
    p::Array{RealVec}
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   PARAMETER STRUCT
#-----------------------------------------------------------------------
"""
    Parameters

Datatype containing all simulation parameters for the post-Minkowskian
N-body integration.

# Fields
- `d::Int`:     Number of spatial dimensions (typically 3)
- `δ::Real`:    Initial time step size
- `rkl::Bool`:  Use the RK4 integrator? If false, uses OrdinaryDiffEq.jl
- `courant::Real`: Courant number for RK4 (adaptive step size)
- `Nrec::Int`:  Skipped timesteps (Nrec=-N keeps last N steps)
- `integrator::Any`: Integrator method for OrdinaryDiffEq.jl
- `atol::Real`: Absolute tolerance for integrators in OrdinaryDiffEq.jl
- `rtol::Real`: Relative tolerance for integrators in OrdinaryDiffEq.jl
- `tspan::Any`: Time span (initial_time, final_time)
- `iter::Int`:  Total number or maximum number of iterations

"""
mutable struct Parameters
    d::Int             # Number of spatial dimensions (typically 3)
    δ::Real            # Initial time step size
    rkl::Bool          # Use rk4? If false, default to OrdinaryDiffEq.jl
    courant::Real      # Courant number for rk4 integrator
    Nrec::Int          # Skipped timesteps (Nrec=-N keeps last N steps)
    integrator::Any    # Integrator method for OrdinaryDiffEq.jl
    atol::Real         # Absolute tolerance in OrdinaryDiffEq.jl
    rtol::Real         # Relative tolerance in OrdinaryDiffEq.jl
    tspan::Any         # Time span (initial_time,final_time)
    iter::Int          # Total number or maximum number of iterations
end #-------------------------------------------------------------------

# Constructor functions for Parameters struct

"""
    Parameters(; kwargs...)

Flexible constructor for Parameters struct with keyword arguments and defaults.

# Keyword Arguments
- `d::Int=3`: Number of spatial dimensions
- `δ::Real=0.01`: Initial time step size for RK4 integrator
- `rkl::Bool=false`: Use RK4 integrator? If false, use OrdinaryDiffEq.jl
- `courant::Real=0.001`: Courant number for RK4 integrator
- `Nrec::Int=100`: Number of timesteps between records
- `integrator::String="Tsit5"`: OrdinaryDiffEq.jl integrator method
- `atol::Real=1e-10`: Absolute tolerance for integrators
- `rtol::Real=1e-10`: Relative tolerance for integrators
- `tspan::Tuple{Real,Real}=(0.0, 10.0)`: Time span for integration
- `iter::Int=10000`: Maximum number of iterations

# Examples
```julia
# Minimal initialization for RK4
params_rk4 = Parameters(rkl=true, tspan=(0.0, 100.0))

# Minimal initialization for Julia integrators
params_julia = Parameters(tspan=(0.0, 50.0), atol=1e-12)

# Full customization
params_custom = Parameters(d=2, δ=0.001, rkl=true, courant=0.0005, 
                          tspan=(0.0, 1000.0), iter=50000)
```
"""
function Parameters(; d::Int=3, δ::Real=0.01, rkl::Bool=false, 
                   courant::Real=0.001, Nrec::Int=100, 
                   integrator::String="Vern8", atol::Real=1e-14, rtol::Real=1e-14,
                   tspan::Tuple{Real,Real}=(0.0, 10.0), iter::Int=10000)
    
    return Parameters(d, δ, rkl, courant, Nrec, integrator, atol, rtol, tspan, iter)
end

"""
    Parameters(tspan::Tuple{Real,Real}; kwargs...)

Convenience constructor that requires only the time span.

# Arguments
- `tspan::Tuple{Real,Real}`: Time span for integration (required)
- `kwargs...`: Additional keyword arguments (see main constructor)

# Example
```julia
params = Parameters((0.0, 100.0), rkl=true, δ=0.005)
```
"""
function Parameters(tspan::Tuple{Real,Real}; kwargs...)
    return Parameters(; tspan=tspan, kwargs...)
end

"""
    ParametersRK4(tspan::Tuple{Real,Real}; kwargs...)

Convenience constructor for RK4 integrator with sensible defaults.

# Arguments
- `tspan::Tuple{Real,Real}`: Time span for integration (required)
- `kwargs...`: Additional keyword arguments

# Example
```julia
params = ParametersRK4((0.0, 1000.0), δ=0.001, courant=0.1)
```
"""
function ParametersRK4(tspan::Tuple{Real,Real}; δ::Real=0.01, 
                       courant::Real=0.1, kwargs...)
    return Parameters(; tspan=tspan, rkl=true, δ=δ, courant=courant, kwargs...)
end

"""
    ParametersJulia(tspan::Tuple{Real,Real}; kwargs...)

Convenience constructor for Julia OrdinaryDiffEq.jl integrators with sensible defaults.

# Arguments
- `tspan::Tuple{Real,Real}`: Time span for integration (required)
- `kwargs...`: Additional keyword arguments

# Example
```julia
params = ParametersJulia((0.0, 50.0), integrator="Vern7", atol=1e-12)
```
"""
function ParametersJulia(tspan::Tuple{Real,Real}; atol::Real=1e-14, rtol::Real=1e-14,
                         integrator::String="Vern8", kwargs...)
    return Parameters(; tspan=tspan, rkl=false, atol=atol, rtol=rtol, integrator=integrator, kwargs...)
end

#-----------------------------------------------------------------------
#   SOLUTION STRUCT
#-----------------------------------------------------------------------
"""
    soln

Mutable datatype for storing integration solutions (only used when RK4 integrator is selected).

# Fields
- `d::Int`: Number of spatial dimensions (typically 3 for 3D space)
- `N::Int`: Number of particles in the system
- `t::RealVec`: Vector recording times at each timestep
- `z::Array{RealVec,1}`: Vector recording phase space coordinates at each timestep
- `zaux::Array{RealVec,1}`: Auxiliary array for intermediate calculations

"""
mutable struct soln
    d::Int                         # number of dimensions
    N::Int                         # number of particles
    t::RealVec                     # Vector recording times at each timestep
    z::Array{RealVec,1}            # Vector recording phase space coordinates at each timestep
    zaux::Array{RealVec,1}
end #-------------------------------------------------------------------
