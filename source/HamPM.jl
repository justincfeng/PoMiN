#-----------------------------------------------------------------------
#   POST-MINKOWSKIAN HAMILTONIAN
#-----------------------------------------------------------------------

module HPM # This module contains the post-Minkowskian Hamiltonian

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

# These functions compute scalar quantities that make up the Hamiltonian

function psf( p::RealVec )
    return dot(p,p)
end

function Enf( m::Real, ps::Real )
    return sqrt(m^2+ps)
end

function rf( qa::RealVec , qb::RealVec )
    return sqrt(dot(qa-qb,qa-qb))
end

function nabf( qa::RealVec , qb::RealVec )
    return (qa-qb)/rf(qa,qb)
end

function ybaf( mb::Real, qb::RealVec , qa::RealVec , pb::RealVec , pa::RealVec )
    return sqrt(mb^2+dot(nabf(qb,qa),pb)^2)/Enf(mb,psf(pb))
end

function Θabf( qa::RealVec , qb::RealVec , pa::RealVec )
    return dot(pa,nabf(qb,qa))
end

function Ξabf( pa::RealVec, pb::RealVec )
    return dot(pa,pb)
end

# Hamiltonian function
function H( d::Int , m::RealVec , Z::RealVec )
    tpfl = typeof(Z[1])
    n = length(m)

    qa, qb, pa, pb  = [zeros(tpfl,d) for _ = 1:4]

    psa, psb, Ena, Enb  = [zero(tpfl) for _ = 1:4]

    rab, yba, Θab, Θba, Ξab  = [zero(tpfl) for _ = 1:5]

    H0, H1, H2, H3  = [zero(tpfl) for _ = 1:4]

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
                yba = ybaf(m[b],qb,qa,pb,pa)
                Θab = Θabf(qa,qb,pa)
                Θba = Θabf(qb,qa,pb)
                Ξab = Ξabf(pa,pb)
                
                H1 -= Ena*Enb*( 1 + psa/(Ena^2) + psb/(Enb^2) )/(2*rab)
                H2 += ( 7*Ξab - Θab*Θba )/(4*rab)

                if m[b] != zero(tpfl)
                    H3 -= 
                    (
                    2*(
                       -2*(Ξab*Θba)^2 + 2*Θab*Θba*Ξab*psb + (Θab*psb)^2 - psb*Ξab^2
                      )/Enb^2
                    +
                    2*(
                       psa*Θba^2 - (Θab*Θba)^2 - 2*Θab*Θba*Ξab + Ξab^2 - psb*Θab^2 
                      )
                    +
                    yba*(
                         3*psa*Θba^2 - (Θab*Θba)^2 - 8*Θab*Θba*Ξab + psa*psb - 3*psb*Θab^2
                        )
                    ) / (4*Ena*Enb*rab*yba*(yba+1)^2)
                # else  # This does not work
                #    H3 -= 
                #    (
                #    4*yba*Ξab^2 - 2*yba*psa*psb + 2*yba*psb*Θab^2 - 3*psa*psb*yba^2 +
                #    yba*psb*Θab^2 - 8*Θab*Θba*Ξab + psa*psb - 3*psb*Θ^2
                #    ) / (4*Ena*Enb*rab*(yba+1)^2)
                end
            end
        end
    end

end


# Derivative of H3 for massless particles  (CHECK THAT THIS IS CAPTURING THE CORRECT CROSS-TERM)
function dH3m0( d::Int , m::RealVec , Z::RealVec )
    tpfl = typeof(Z[1])
    n = length(m)
    nn = 2*n*d

    o = zero(tpfl)

    qa, qb, qc, pa, pb, pc  = [zeros(tpfl,d) for _ = 1:6]

    psa, psb, psc, Ena, Enb, Enc = [o for _ = 1:6]

    rcb, ybc, Θcb, Θbc, Ξcb  = [o for _ = 1:5]
    rac, yca, Θac, Θca, Ξac  = [o for _ = 1:5]

    dpsa, dEna, dEna = [o for _ = 1:3]
    dpsb, dEnb, dEnb = [o for _ = 1:3]
    dpsc, dEnc, dEnc = [o for _ = 1:3]

    drcb, dybc, dΘcb, dΘbc, dΞcb  = [o for _ = 1:5]
    drac, dyca, dΘac, dΘca, dΞac  = [o for _ = 1:5]

    dH3 = zeros(tpfl,nn)

    for c=1:n

    if m[c] == zero(tpfl)

        qc = Z2q(n,d,c,Z)
        pc = Z2p(n,d,c,Z)
        psc = psf(pc)
        Enc = Enf(m[c],psb)

        for i=1:d
            for a=1:n

                qa = Z2q(n,d,a,Z)
                pa = Z2p(n,d,a,Z)

                psa = psf(pa)
                Ena = Enf(m[a],psa)

                rac = rf(qa,qc)
                yca = ybaf(m[c],qc,qa,pc,pa)
                Θac = Θabf(qa,qc,pa)
                Θca = Θabf(qc,qa,pc)
                Ξac = Ξabf(pa,pc)

                # qs
                dpsa = o
                dpsc = o
                dEna = o
                dEnc = o
                drac = (qc[i]-qa[i])/rac
                dΘac = (pc[i]-Θac*drac)/rac
                dΞac = o

                dyac = sign(Θac)*(Ena*dΘac - dEna*Θac)/(Ena^2)

                dH3[d*(c-1)+i] -= 
                        (1/4)*(-((2*Ena*Enc*dyca*rac*yca*(2*psc*psc*Θac*Θac + 4*psc*Θac*Θca*Ξac - 2*(psc - 2*Θca*Θca)*Ξac*Ξac + Enc*Enc*(2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(Θac*Θca - 8*Ξac)*yca - psa*Θca*Θca*(2 + 3*yca) + psc*(psa*yca - Θac*Θac*(2 + 3*yca)))) + dEnc*Ena*rac*yca*(1 + yca)*(2*psc*psc*Θac*Θac + 4*psc*Θac*Θca*Ξac - 2*(psc - 2*Θca*Θca)*Ξac*Ξac + Enc*Enc*(2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(Θac*Θca - 8*Ξac)*yca - psa*Θca*Θca*(2 + 3*yca) + psc*(psa*yca - Θac*Θac*(2 + 3*yca)))) - Ena*Enc*dyca*rac*(1 + yca)*(-2*psc*psc*Θac*Θac - 4*Θca*Θca*Ξac*Ξac + 2*psc*Ξac*(-2*Θac*Θca + Ξac) + Enc*Enc*(-2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(-(Θac*Θca) + 8*Ξac)*yca + psa*Θca*Θca*(2 + 3*yca) + psc*(-(psa*yca) + Θac*Θac*(2 + 3*yca)))) - Ena*Enc*drac*yca*(1 + yca)*(-2*psc*psc*Θac*Θac - 4*Θca*Θca*Ξac*Ξac + 2*psc*Ξac*(-2*Θac*Θca + Ξac) + Enc*Enc*(-2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(-(Θac*Θca) + 8*Ξac)*yca + psa*Θca*Θca*(2 + 3*yca) + psc*(-(psa*yca) + Θac*Θac*(2 + 3*yca)))) - dEna*Enc*rac*yca*(1 + yca)*(-2*psc*psc*Θac*Θac - 4*Θca*Θca*Ξac*Ξac + 2*psc*Ξac*(-2*Θac*Θca + Ξac) + Enc*Enc*(-2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(-(Θac*Θca) + 8*Ξac)*yca + psa*Θca*Θca*(2 + 3*yca) + psc*(-(psa*yca) + Θac*Θac*(2 + 3*yca)))) + Ena*rac*yca*(1 + yca)*(2*dpsc*Enc*Enc*Enc*Θac*Θac + 4*psc*psc*Θac*(-(Enc*dΘac) + dEnc*Θac) - 4*Enc*Enc*Enc*dΞac*Ξac + 4*Enc*Enc*Enc*dΘca*Θac*Ξac + 2*dpsc*Enc*Ξac*Ξac - dpsc*Enc*Enc*Enc*psa*yca + 3*dpsc*Enc*Enc*Enc*Θac*Θac*yca + 8*Enc*Enc*Enc*dΘca*Θac*Ξac*yca + Θca*Θca*(-8*Enc*dΞac*Ξac + 8*dEnc*Ξac*Ξac + Enc*Enc*Enc*(dyca*(3*psa - Θac*Θac) - 2*dΘac*Θac*(2 + yca) + dpsa*(2 + 3*yca))) + psc*((-4*dpsc*Enc + 3*Enc*Enc*Enc*dyca)*Θac*Θac + 4*Enc*(dΞac - dΘac*Θca)*Ξac - 4*dEnc*Ξac*Ξac - Enc*Enc*Enc*(psa*dyca + dpsa*yca) + Θac*(-4*Enc*dΘca*Ξac + Θca*(-4*Enc*dΞac + 8*dEnc*Ξac) + 2*Enc*Enc*Enc*dΘac*(2 + 3*yca))) + Θca*(-4*Enc*Ξac*(dpsc*Θac + 2*dΘca*Ξac) + Enc*Enc*Enc*(4*dΞac*Θac*(1 + 2*yca) + 4*Ξac*(2*dyca*Θac + dΘac*(1 + 2*yca)) + 2*dΘca*(-(Θac*Θac*(2 + yca)) + psa*(2 + 3*yca))))))/(Ena*Ena*Enc*Enc*Enc*Enc*rac*rac*yca*yca*(1 + yca)*(1 + yca)*(1 + yca))))

                # ps
                dpsc = 2*pc
                dEnc = dpsc/(2*Enc)
                if c==a
                    dpsa = dpsc
                else
                    dpsa = o
                end
                dEna = dpsa/(2*Ena)
                drac = o
                dΘac = o
                dΞac = pa[i]

                dyac = sign(Θac)*(Ena*dΘac - dEna*Θac)/(Ena^2)

                dH3[d*(c-1+n)+i] -= 
                        (1/4)*(-((2*Ena*Enc*dyca*rac*yca*(2*psc*psc*Θac*Θac + 4*psc*Θac*Θca*Ξac - 2*(psc - 2*Θca*Θca)*Ξac*Ξac + Enc*Enc*(2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(Θac*Θca - 8*Ξac)*yca - psa*Θca*Θca*(2 + 3*yca) + psc*(psa*yca - Θac*Θac*(2 + 3*yca)))) + dEnc*Ena*rac*yca*(1 + yca)*(2*psc*psc*Θac*Θac + 4*psc*Θac*Θca*Ξac - 2*(psc - 2*Θca*Θca)*Ξac*Ξac + Enc*Enc*(2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(Θac*Θca - 8*Ξac)*yca - psa*Θca*Θca*(2 + 3*yca) + psc*(psa*yca - Θac*Θac*(2 + 3*yca)))) - Ena*Enc*dyca*rac*(1 + yca)*(-2*psc*psc*Θac*Θac - 4*Θca*Θca*Ξac*Ξac + 2*psc*Ξac*(-2*Θac*Θca + Ξac) + Enc*Enc*(-2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(-(Θac*Θca) + 8*Ξac)*yca + psa*Θca*Θca*(2 + 3*yca) + psc*(-(psa*yca) + Θac*Θac*(2 + 3*yca)))) - Ena*Enc*drac*yca*(1 + yca)*(-2*psc*psc*Θac*Θac - 4*Θca*Θca*Ξac*Ξac + 2*psc*Ξac*(-2*Θac*Θca + Ξac) + Enc*Enc*(-2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(-(Θac*Θca) + 8*Ξac)*yca + psa*Θca*Θca*(2 + 3*yca) + psc*(-(psa*yca) + Θac*Θac*(2 + 3*yca)))) - dEna*Enc*rac*yca*(1 + yca)*(-2*psc*psc*Θac*Θac - 4*Θca*Θca*Ξac*Ξac + 2*psc*Ξac*(-2*Θac*Θca + Ξac) + Enc*Enc*(-2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + Θac*Θca*(-(Θac*Θca) + 8*Ξac)*yca + psa*Θca*Θca*(2 + 3*yca) + psc*(-(psa*yca) + Θac*Θac*(2 + 3*yca)))) + Ena*rac*yca*(1 + yca)*(2*dpsc*Enc*Enc*Enc*Θac*Θac + 4*psc*psc*Θac*(-(Enc*dΘac) + dEnc*Θac) - 4*Enc*Enc*Enc*dΞac*Ξac + 4*Enc*Enc*Enc*dΘca*Θac*Ξac + 2*dpsc*Enc*Ξac*Ξac - dpsc*Enc*Enc*Enc*psa*yca + 3*dpsc*Enc*Enc*Enc*Θac*Θac*yca + 8*Enc*Enc*Enc*dΘca*Θac*Ξac*yca + Θca*Θca*(-8*Enc*dΞac*Ξac + 8*dEnc*Ξac*Ξac + Enc*Enc*Enc*(dyca*(3*psa - Θac*Θac) - 2*dΘac*Θac*(2 + yca) + dpsa*(2 + 3*yca))) + psc*((-4*dpsc*Enc + 3*Enc*Enc*Enc*dyca)*Θac*Θac + 4*Enc*(dΞac - dΘac*Θca)*Ξac - 4*dEnc*Ξac*Ξac - Enc*Enc*Enc*(psa*dyca + dpsa*yca) + Θac*(-4*Enc*dΘca*Ξac + Θca*(-4*Enc*dΞac + 8*dEnc*Ξac) + 2*Enc*Enc*Enc*dΘac*(2 + 3*yca))) + Θca*(-4*Enc*Ξac*(dpsc*Θac + 2*dΘca*Ξac) + Enc*Enc*Enc*(4*dΞac*Θac*(1 + 2*yca) + 4*Ξac*(2*dyca*Θac + dΘac*(1 + 2*yca)) + 2*dΘca*(-(Θac*Θac*(2 + yca)) + psa*(2 + 3*yca))))))/(Ena*Ena*Enc*Enc*Enc*Enc*rac*rac*yca*yca*(1 + yca)*(1 + yca)*(1 + yca))))
            end
            for b=1:n

                qb = Z2q(n,d,b,Z)
                pb = Z2p(n,d,b,Z)

                psb = psf(pb)
                Enb = Enf(m[b],psb)

                rcb = rf(qc,qb)
                ybc = ybaf(m[b],qb,qc,pb,pc)
                Θcb = Θabf(qc,qb,pc)
                Θbc = Θabf(qb,qc,pb)
                Ξcb = Ξabf(pc,pb)

                # qs
                dpsb = o
                dpsc = o
                dEnb = o
                dEnc = o
                drcb = -(qb[i]-qc[i])/rcb
                dΘcb = (pc[i]-Θbc*drcb)/rcb
                dΘbc = (pb[i]+Θcb*drcb)/rcb
                dΞcb = o

                dycb = sign(Θcb)*(Enc*dΘcb - dEnc*Θcb)/(Enc^2)

                signTh = sign(Θbc)

                dH3[d*(c-1)+i] -= 
                        (1/4)*(2*Enc*psb*drcb*(1 + ybc)*(8*signTh*sqrt(psb)*Θcb*Ξcb*ybc - 4*Ξcb*Ξcb*ybc - psb*Θcb*Θcb*(-1 + ybc)*(3 + ybc) +  psc*psb*(1 + ybc)*(-1 + 3*ybc)) + rcb*(dpsb*Enc*(1 + ybc)*(8*signTh*sqrt(psb)*Θcb*Ξcb*ybc - 2*6*Ξcb*Ξcb*ybc - psb*Θcb*Θcb*(3 + ybc*(2 + ybc)) + psc*psb*(1 + ybc*(2 + 3*ybc))) - 2*(dpsc*Enc*psb*psb*(1 + ybc)*(1 + ybc)*(-1 + 3*ybc) + dEnc*psb*(1 + ybc)*(-8*signTh*sqrt(psb)*Θcb*Ξcb*ybc + 4*Ξcb*Ξcb*ybc + psb*Θcb*Θcb*(-1 + ybc)*(3 + ybc) - psc*psb*(1 + ybc)*(-1 + 3*ybc)) + 2*Enc*sqrt(psb)*(dΘbc*(1 + ybc)*(4*sqrt(psb)*Θcb*Ξcb - 4*signTh*Ξcb*Ξcb - signTh*psb*Θcb*Θcb*(2 + ybc) + signTh*psc*psb*(2 + 3*ybc)) + sqrt(psb)*(4*dΞcb*(signTh*sqrt(psb)*Θcb - Ξcb)*ybc*(1 + ybc) - sqrt(psb)*dΘcb*(1 + ybc)*(-4*signTh*Ξcb*ybc + sqrt(psb)*Θcb*(-1 + ybc)*(3 + ybc)) + dybc*(-8*signTh*sqrt(psb)*Θcb*Ξcb*ybc + 2*Ξcb*Ξcb*(1 + 3*ybc) + psb*(-3*psc*ybc*(1 + ybc) + Θcb*Θcb*(-2 + ybc*(3 + ybc)))))))))/(2*Enc*Enc*sqrt(psb)*psb*rcb*rcb*(1 + ybc)*(1 + ybc)*(1 + ybc))
                # ps
                dpsc = 2*pc
                dEnc = dpsc/(2*Enc)
                if c==b
                    dpsb = dpsc
                else
                    dpsb = o
                end
                dEnb = dpsb/(2*Enb)
                drbc = o
                dΘcb = (qb[i]-qc[i])/rcb
                dΘbc = (qc[i]-qb[i])/rcb
                dΞcb = pb[i]

                dycb = sign(Θcb)*(Enc*dΘcb - dEnc*Θcb)/(Enc^2)

                signTh = sign(Θbc)

                dH3[d*(c-1+n)+i] -= 
                        (1/4)*(2*Enc*psb*drcb*(1 + ybc)*(8*signTh*sqrt(psb)*Θcb*Ξcb*ybc - 4*Ξcb*Ξcb*ybc - psb*Θcb*Θcb*(-1 + ybc)*(3 + ybc) +  psc*psb*(1 + ybc)*(-1 + 3*ybc)) + rcb*(dpsb*Enc*(1 + ybc)*(8*signTh*sqrt(psb)*Θcb*Ξcb*ybc - 2*6*Ξcb*Ξcb*ybc - psb*Θcb*Θcb*(3 + ybc*(2 + ybc)) + psc*psb*(1 + ybc*(2 + 3*ybc))) - 2*(dpsc*Enc*psb*psb*(1 + ybc)*(1 + ybc)*(-1 + 3*ybc) + dEnc*psb*(1 + ybc)*(-8*signTh*sqrt(psb)*Θcb*Ξcb*ybc + 4*Ξcb*Ξcb*ybc + psb*Θcb*Θcb*(-1 + ybc)*(3 + ybc) - psc*psb*(1 + ybc)*(-1 + 3*ybc)) + 2*Enc*sqrt(psb)*(dΘbc*(1 + ybc)*(4*sqrt(psb)*Θcb*Ξcb - 4*signTh*Ξcb*Ξcb - signTh*psb*Θcb*Θcb*(2 + ybc) + signTh*psc*psb*(2 + 3*ybc)) + sqrt(psb)*(4*dΞcb*(signTh*sqrt(psb)*Θcb - Ξcb)*ybc*(1 + ybc) - sqrt(psb)*dΘcb*(1 + ybc)*(-4*signTh*Ξcb*ybc + sqrt(psb)*Θcb*(-1 + ybc)*(3 + ybc)) + dybc*(-8*signTh*sqrt(psb)*Θcb*Ξcb*ybc + 2*Ξcb*Ξcb*(1 + 3*ybc) + psb*(-3*psc*ybc*(1 + ybc) + Θcb*Θcb*(-2 + ybc*(3 + ybc)))))))))/(2*Enc*Enc*sqrt(psb)*psb*rcb*rcb*(1 + ybc)*(1 + ybc)*(1 + ybc))
            end
        end

    end   # End massless check
    end   # End c loop

    return dH3

end

function dH( d::Int , m::RealVec , Z::RealVec )
    return ForwardDiff.gradient(x->H(d,m,x),Z) + dH3m0( d , m , Z )
end

end # module HPM