const RealVec{T<:Real} = Array{T,1}

struct Particles
    m::RealVec
    q::Array{RealVec,1}
    p::Array{RealVec,1}
end

struct Parameters
    sym :: Bool
    rkl :: Tuple{Bool,Real}
    jli :: Bool
    ω :: Real
    tspan :: Tuple{Real,Real}
    iter  :: Int
    tol   :: Real
end


struct soln
    t::RealVec
    z::Array{RealVec,1}
end