module HamHO

using LinearAlgebra
using ForwardDiff

include("pomin-types.jl")

function H( z::RealVec )
    tpfl=typeof(z[1])
    n2 = length(z)

    h = zero(tpfl)

    if iseven(n2)
        n = Int(n2 / 2)

        for i=1:n
            q = z[i]
            p = z[n+i]
            h += dot(p,p)/2 + dot(q,q)/2
        end
    
        return h
    else
        return h
    end
end

function dH( z::RealVec )
    return ForwardDiff.gradient(H,z)
end

end # end of Hamiltonian

