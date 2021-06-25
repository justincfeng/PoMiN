module exf      # This module contains external forces

using LinearAlgebra
using ForwardDiff

const RealVec{T<:Real} = Array{T,1}

# These functions extract particle positions and momenta from Z

function Z2q( n::Int , d::Int , i::Int , Z::RealVec )
    tpfl = typeof(Z[1])
    q = zeros(tpfl,d)
    for j=1:d
        q[j] = Z[d*(i-1)+j]
    end
    return q
end

function Z2p( n::Int , d::Int , i::Int , Z::RealVec )
    tpfl = typeof(Z[1])
    p = zeros(tpfl,d)
    for j=1:d
        p[j] = Z[d*(i-1+n)+j]
    end
    return p
end

# Hamiltonian function
function H( d::Int , m::RealVec , Z::RealVec )
    tpfl = typeof(Z[1])
    n = length(m)

    qa, pa = [zeros(tpfl,d) for _ = 1:2]

    Hx  = zero(tpfl)

    for a=1:n
        qa = Z2q(n,d,a,Z)
        pa = Z2p(n,d,a,Z)

        Hx += zero(tpfl)
    end

end

function dH( d::Int , m::RealVec , Z::RealVec , H::Function )
    return ForwardDiff.gradient(H,Z)
end

end # module HPM