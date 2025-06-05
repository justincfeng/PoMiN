

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