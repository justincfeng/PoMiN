#-----------------------------------------------------------------------
#   GRAVITATIONAL WAVE STRAIN CALCULATOR
#-----------------------------------------------------------------------
module gwsc      # This module contains functions for computing gw strain

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

# This function computes the linearized gw strain for a system of particles
function hstr( d::Int , m::RealVec , Z::RealVec , θ::Real , ϕ::Real )
    tpfl = typeof(Z[1])
    n = length(m)

    q, p = [zeros(tpfl,d) for _ = 1:2]

    hstr = zero(tpfl)
    En = zero(tpfl)

    St = sin(θ)
    Sp = sin(ϕ)
    Ct = cos(θ)
    Cp = cos(ϕ)

    for a=1:n
        q = Z2q(n,d,a,Z)
        p = Z2p(n,d,a,Z)

        En = sqrt(m[a]^2 + dot(p,p))

        hp += 
            (
                -8*Ct*p[2]*p[3]*Sp*St + 2*(p[3]^2)*(1 - (Ct^2) + (St^2)) + 
                (p[1]^2)*(
                            -1 - 3*(Sp^2) - (Ct^2)*(-1 + (Sp^2)) + 
                            (-1 + (Sp^2))*(St^2) + (Cp^2)*(3 + (Ct^2) - (St^2))
                         ) +
                4*Cp*p[1]*(-2*Ct*p[3]*St + p[2]*Sp*(3 + (Ct^2) - (St^2))) + 
                (p[2]^2)*(
                            -1 + 3*(Sp^2) + (Ct^2)*(1 + (Sp^2)) - 
                            (1 + (Sp^2))*(St^2) + (Cp^2)*(-3 - (Ct^2) + (St^2))
                         )
            )/(2*En)

        hx += (4*(Cp*p[2] - p[1]*Sp)*(Cp*Ct*p[1] + Ct*p[2]*Sp - p[3]*St))/En
    end

    return [ hp ; hx ]
end

end # module HPM