"""
    RealVec{T<:Real}

Type alias for a vector of real numbers. Equivalent to `Array{T}` where `T` is a subtype of `Real`.
Used throughout PoMiN for representing position and momentum vectors in phase space.
"""
const RealVec{T<:Real} = Array{T}

"""
    RealMtx{T<:Real}

Type alias for a matrix of real numbers. Equivalent to `Array{T,2}` where `T` is a subtype of `Real`.
Used for representing transformation matrices and other 2D arrays in PoMiN calculations.
"""
const RealMtx{T<:Real} = Array{T,2}

"""
    Particles

Datatype representing a collection of particles in the N-body system.

# Fields
- `m::RealVec`: Vector of particle masses
- `q::Array{RealVec}`: Array of position vectors for each particle
- `p::Array{RealVec}`: Array of momentum vectors for each particle

# Notes
This structure is used to store the physical properties and state of all particles
in the post-Minkowskian N-body simulation. Each particle has a mass, position vector,
and momentum vector in the relativistic framework.
"""
struct Particles
    m::RealVec
    q::Array{RealVec}
    p::Array{RealVec}
end #-------------------------------------------------------------------

"""
    Parameters

Datatype containing all simulation parameters for the post-Minkowskian N-body integration.

# Fields
- `d::Int`: Number of spatial dimensions (typically 3 for 3D space)
- `δ::Real`: Initial time step size for RK4 integrator
- `rkl::Tuple{Bool}`: Use the RK4 integrator? If false, defaults to OrdinaryDiffEq.jl integrator
- `courant::Real`: Courant number for RK4 integrator (controls timestep stability)
- `Nrec::Int`: Number of timesteps between records (Nrec=-N records the last N timesteps)
- `integrator::String`: Integrator method for OrdinaryDiffEq.jl (e.g., "Vern7", "Rodas5")
- `atol::Real`: Absolute tolerance for integrators in OrdinaryDiffEq.jl (controls accuracy)
- `rtol::Real`: Relative tolerance for integrators in OrdinaryDiffEq.jl (controls accuracy)
- `tspan::Tuple{Real,Real}`: Time span as (initial_time, final_time)
- `iter::Int`: Number of iterations or maximum iterations, depending on integrator

# Notes
This structure encapsulates all the numerical parameters needed to configure
the integration of Hamilton's equations in the post-Minkowskian approximation.
The choice between RK4 and OrdinaryDiffEq.jl integrators allows for different
balances between speed and accuracy.
"""
mutable struct Parameters
    d::Int                         # Number of spatial dimensions (typically 3)
    δ::Real                        # Initial time step size
    rkl::Bool                      # Use the rk4 integrator? If false, default to OrdinaryDiffEq.jl integrator
    courant::Real                  # Courant number for rk4 integrator
    Nrec::Int                      # Number of timesteps between records (Nrec=-N records the last N timesteps)
    integrator::String             # Integrator method for OrdinaryDiffEq.jl
    atol::Real                     # Absolute tolerance for integrators in OrdinaryDiffEq.jl
    rtol::Real                     # Relative tolerance for integrators in OrdinaryDiffEq.jl
    tspan::Tuple{Real,Real}        # Time span (initial_time,final_time)
    iter::Int                      # Either the number of iterations, or maximum number of iterations, depending on integrator
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
                   integrator::String="Tsit5", atol::Real=1e-10, rtol::Real=1e-10,
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
params = ParametersRK4((0.0, 1000.0), δ=0.001, courant=0.0005)
```
"""
function ParametersRK4(tspan::Tuple{Real,Real}; δ::Real=0.01, 
                       courant::Real=0.001, kwargs...)
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
function ParametersJulia(tspan::Tuple{Real,Real}; atol::Real=1e-10, rtol::Real=1e-10,
                         integrator::String="Tsit5", kwargs...)
    return Parameters(; tspan=tspan, rkl=false, atol=atol, rtol=rtol, integrator=integrator, kwargs...)
end

"""
    soln

Mutable datatype for storing integration solutions (only used when RK4 integrator is selected).

# Fields
- `d::Int`: Number of spatial dimensions (typically 3 for 3D space)
- `N::Int`: Number of particles in the system
- `t::RealVec`: Vector recording times at each timestep
- `z::Array{RealVec,1}`: Vector recording phase space coordinates at each timestep
- `zaux::Array{RealVec,1}`: Auxiliary array for intermediate calculations

# Notes
This mutable structure is used to store the time evolution of the N-body system
when using the custom RK4 integrator. The phase space coordinates `z` contain
both positions and momenta for all particles at each recorded timestep.
The structure is mutable to allow efficient updating during integration.
"""
mutable struct soln
    d::Int                         # number of dimensions
    N::Int                         # number of particles
    t::RealVec                     # Vector recording times at each timestep
    z::Array{RealVec,1}            # Vector recording phase space coordinates at each timestep
    zaux::Array{RealVec,1}
end #-------------------------------------------------------------------

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

# Notes
This utility function is used throughout PoMiN to generate type-consistent
numeric constants, enabling support for different precision arithmetic
(Float64, BigFloat, DoubleFloats, etc.) in post-Minkowskian calculations.
If the input type is not a subtype of Real, it defaults to Float64.
"""
function tpnum(tpfl::Type)
    if tpfl <: Real
        return (tpfl(0),tpfl(1),tpfl(2),tpfl(3),tpfl(4),
                tpfl(5),tpfl(6),tpfl(7),tpfl(8),tpfl(9))
    else
        return (Float64(0.0),Float64(1.0),Float64(2.0),Float64(3.0),
                Float64(4.0),Float64(5.0),Float64(6.0),Float64(7.0),
                Float64(8.0),Float64(9.0))
    end
end #-------------------------------------------------------------------