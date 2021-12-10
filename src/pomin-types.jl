const RealVec{T<:Real} = Array{T}

#   Particle datatype
struct Particles
    m::RealVec
    q::Array{RealVec}
    p::Array{RealVec}
end

#   Datatype for parameters
struct Parameters
    sym :: Tuple{Bool,Real}     # Use the symplectic integrator?            The real parameter in the tuple is the ω parameter in the Tao map
    rkl :: Tuple{Bool,Real}     # Use the rk4 integrator?                   The real parameter in the tuple is the Courant number
    jli :: Tuple{Bool,Real}     # Use integrators in OrdinaryDiffEq.jl?     The real parameter in the tuple is the tolerance
    tspan :: Tuple{Real,Real}   # Time span (initial_time,final_time)
    iter  :: Int                # Either the number of iterations, or maximum number of iterations, depending on integrator
end

#   Datatype for solutions
struct soln
    t::RealVec              # Vector recording times at each timestep
    z::Array{RealVec,1}     # Vector recording phase space coordinates at each timestep
    zaux::Array{RealVec,1}
end