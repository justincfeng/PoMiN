const RealVec{T<:Real} = Array{T}
const RealMtx{T<:Real} = Array{T,2}



#   Particle datatype
struct Particles
    m::RealVec
    q::Array{RealVec}
    p::Array{RealVec}
end #-------------------------------------------------------------------

#   Datatype for parameters
struct Parameters
    sym::Tuple{Bool,Real}     # Use the symplectic integrator?            The real parameter in the tuple is the ω parameter in the Tao map
    rkl::Tuple{Bool,Real}     # Use the rk4 integrator?                   The real parameter in the tuple is the Courant number
    jli::Tuple{Bool,Real}     # Use integrators in OrdinaryDiffEq.jl?     The real parameter in the tuple is the tolerance
    tspan::Tuple{Real,Real}   # Time span (initial_time,final_time)
    iter::Int                # Either the number of iterations, or maximum number of iterations, depending on integrator
end #-------------------------------------------------------------------

#   Datatype for solutions
mutable struct soln
    d::Int                  # number of dimensions
    N::Int                  # number of particles
    t::RealVec              # Vector recording times at each timestep
    z::Array{RealVec,1}     # Vector recording phase space coordinates at each timestep
    zaux::Array{RealVec,1}
end #-------------------------------------------------------------------

#   Generate constants for a given datatype
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

#-----------------------------------------------------------------------
"""
    to_phase_tuple(system::Particles)

Convert a Particles object into a tuple of (Z, m, d) suitable for use with FHE.

Arguments:
- system: A Particles object

Returns:
- Z: Phase space vector [q1...qN, p1...pN]
- m: Mass vector
- d: Number of dimensions (3)

Example:
```julia
Z, m, d = to_phase_tuple(system)
FHE(Z, m, d)  # Ready for Hamilton's equations
```
"""
function to_phase_tuple(system::Particles)
    # Get number of dimensions from first position vector
    d = length(system.q[1])
    N = length(system.m)  # number of particles
    tpfl = eltype(system.m)  # system's numeric type
    
    # Initialize phase space vector
    Z = zeros(tpfl, 2*N*d)
    
    # Fill positions (first half)
    for i in 1:N
        for j in 1:d
            Z[d*(i-1) + j] = system.q[i][j]
        end
    end
    
    # Fill momenta (second half)
    for i in 1:N
        for j in 1:d
            Z[d*(i-1+N) + j] = system.p[i][j]
        end
    end
    
    return Z, copy(system.m), d
end #-------------------------------------------------------------------