#-----------------------------------------------------------------------
#   POST-MINKOWSKIAN HAMILTONIAN
#-----------------------------------------------------------------------

module HamPM

using LinearAlgebra
using ForwardDiff

include("../../pomin-types.jl")
include("HamTools.jl")

# These functions compute scalar quantities that make up the Hamiltonian

"""
    psf( p::RealVec )

Momentum squared function. Returns ``p^2=\\vec{p}⋅\\vec{p}`` given
relativistic momentum ``\\vec{p}``.
"""
function psf( p::RealVec )
    return dot(p,p)
end #-------------------------------------------------------------------

"""
    Enf( m::Real, ps::Real )

Kinetic energy function. Returns ``E=\\sqrt{m^2+p^2}`` for massive
particles (m>0), and ``E=\\sqrt{p^2}=|p|`` for massless particles (m=0),
given ``m`` and ``\\vec{p}``.
"""
function Enf( m::Real, ps::Real )
    if m > zero(m)
        return sqrt(m^2+ps)
    else
        return sqrt(ps)
    end
end #-------------------------------------------------------------------

"""
    rf( qa::RealVec , qb::RealVec )

Euclidean distance function. Returns ``r_{ab}=\\sqrt{(qa-qb)⋅(qa-qb)}``,
given particle positions ``\\vec{q}_a`` and ``\\vec{q}_b``.
"""
function rf( qa::RealVec , qb::RealVec )
    return norm(qa-qb)
end #-------------------------------------------------------------------

"""
    nabf( qa::RealVec , qb::RealVec )

Separation unit vector. Returns
``\\vec{n}_{ab}=(\\vec{q}_a-\\vec{q}_b)/\\sqrt{r_{ab}}``, given particle
positions ``\\vec{q}_a`` and ``\\vec{q}_b``.
"""
function nabf( qa::RealVec , qb::RealVec )
    return (qa-qb)/rf(qa,qb)
end #-------------------------------------------------------------------

"""
    ybaf( mb::Real, qb::RealVec , qa::RealVec , pb::RealVec ,)

Returns ``y_{ba}=\\sqrt{m_b^2+(\\vec{n}_{ab}⋅\\vec{p}_b)^2}/E_b``, where
``E_b`` is Kinetic energy for particle ``b``. Inputs are particle ``b``
mass ``m_b``, particle positions ``\\vec{q}_b`` and ``\\vec{q}_a``, and
particle momenta ``\\vec{p}_b`` and ``\\vec{p}_a``.
"""
function ybaf( mb::Real, qb::RealVec , qa::RealVec , pb::RealVec )
    return sqrt(mb^2+Θabf(qb,qa,pb)^2)/Enf(mb,psf(pb))
end #-------------------------------------------------------------------

"""
    Θabf( qa::RealVec , qb::RealVec , pa::RealVec )

Returns ``Θ_{ba}=\\vec{p}_a⋅\\vec{n}_{ab}``, where ``\\vec{p}_b`` is
momentum for particle ``b``, and ``\\vec{n}_{ab}`` is unit separation
vector.
"""
function Θabf( qa::RealVec , qb::RealVec , pa::RealVec )
    return dot(pa,nabf(qb,qa))
end #-------------------------------------------------------------------

"""
    Ξabf( pa::RealVec, pb::RealVec )

Returns ``Ξ_{ba}=\\vec{p}_a⋅\\vec{p}_{b}``, where ``\\vec{p}_b`` and
``\\vec{p}_a`` are momenta for particles ``a`` and ``b``.
"""
function Ξabf( pa::RealVec, pb::RealVec )
    return dot(pa,pb)
end #-------------------------------------------------------------------

"""
    H( Z::RealVec , m::RealVec , d::Int = 3 )

Hamiltonian function. Returns the Hamiltonian for a system of particles
with masses ``m`` and positions ``Z``.
"""
function H( Z::RealVec , m::RealVec , d::Int = 3 )
    tpfl = typeof(Z[1])

    ν0, ν1, ν2, ν3, ν4, ν5, ν6, ν7, ν8, ν9 = tpnum(tpfl)

    n = length(m)

    qa, qb, pa, pb  = [zeros(tpfl,d) for _ = 1:4]

    psa, psb, Ena, Enb  = [ν0 for _ = 1:4]

    rab, yba, Θab, Θba, Ξab  = [ν0 for _ = 1:5]

    H0, H1, H2, H3  = [ν0 for _ = 1:4]

    for a=1:n
        qa = Z2q(n,d,a,Z)
        pa = Z2p(n,d,a,Z)

        psa = psf(pa)
        Ena = Enf(m[a],psa)

        H0 += Ena

        for b=1:n
            if b!=a
                qb = Z2q(n,d,b,Z)
                pb = Z2p(n,d,b,Z)

                psb = psf(pb)
                Enb = Enf(m[b],psb)

                rab = rf(qa,qb)
                yba = ybaf(m[b],qb,qa,pb)
                Θab = Θabf(qa,qb,pa)
                Θba = Θabf(qb,qa,pb)
                Ξab = Ξabf(pa,pb)
                
                H1 -= Ena*Enb*( ν1 + psa/(Ena^2) + psb/(Enb^2) )/(ν2*rab)
                H2 += ( ν7*Ξab - Θab*Θba )/(ν4*rab)

                if m[b] != ν0 && yba != ν0 && yba != -ν1
                    # H3 massive case - corrected to match C implementation
                    H3 -= (ν1/(ν4*rab*Ena*Enb*yba*(yba+ν1)^2)) * (
                        ν2 * ( ν2*Ξab^2*Θba^2 + ν2*Θab*Θba*Ξab*psb + Θab^2*psb^2 - Ξab^2*psb ) /(Enb^2) +
                        ν2 * (-psa*Θba^2 + Θab^2*Θba^2 - ν2*Θab*Θba*Ξab + Ξab^2 - Θab^2*psb) +
                        ((-ν3*psa*Θba^2 + Θab^2*Θba^2 - ν8*Θab*Θba*Ξab + psa*psb - ν3*Θab^2*psb) * yba)
                    )
                else
                    H3 += (ν3*Enb*psb*Θab^2 - Enb*psa*(psb - ν3*Θba^2) + 
                           Enb*Θab*Θba*(-(Θab*Θba) + ν8*Ξab) + ν2*psa*psb*abs(Θba) - 
                           ν2*psb*Θab^2*abs(Θba) - ν4*Ξab^2*abs(Θba)
                           ) / (ν4*Ena*rab*(Enb + abs(Θba))^2)
                end
            end
        end
    end

    return H0+H1+H2+H3

end #-------------------------------------------------------------------

"""
    dH( Z::RealVec , m::RealVec , d::Int = 3 )

Gradient of the Hamiltonian function.
"""
function dH( Z::RealVec , m::RealVec , d::Int = 3 )
    return ForwardDiff.gradient(x -> H(x, m, d), Z)
end #-------------------------------------------------------------------

# Right hand side of Hamilton's equations
"""
    FHE( Z::RealVec , m::RealVec , d::Int = 3 )

Right hand side of Hamilton's equations.
"""
FHE = (Z,m,d=3)->Jsympl(dH(Z,m,d))
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    dp_scatter(p, b, m1, m2; G=1, c=1)

Analytical momentum exchange formula for post-Minkowskian scattering.
Returns the change in momentum for particle 1 in the y-direction.
"""
function dp_scatter(p, b, m1, m2; Ga=1, ca=1)
    tpfl = typeof(p)
    ν0, ν1, ν2, ν3, ν4, ν5, ν6, ν7, ν8, ν9 = tpnum(tpfl)
    c = tpfl(ca)
    G = tpfl(Ga)

    E1 = c*sqrt((c*m1)^2 + p^2)
    E2 = c*sqrt((c*m2)^2 + p^2)
    # Formula matches LSB_momX.txt exactly
    dp = (ν2*G*(E1*E2)^2)/(p*(E1+E2)) *
         (ν1 + (ν1/E1^2 + ν1/E2^2 + ν4/(E1*E2))*p^2 + p^4/(E1*E2)^2) / b
    return dp
end #-------------------------------------------------------------------

export dp_scatter

end # end of HamPM
