#-----------------------------------------------------------------------
#   POST-MINKOWSKIAN HAMILTONIAN
#-----------------------------------------------------------------------

# These functions compute scalar quantities that make up the Hamiltonian

"""
    psf( p::RealVec )

Momentum squared function. Returns ``p^2=\\vec{p}⋅\\vec{p}`` given 
relativistic momentum ``\\vec{p}``.
"""
function psf( p::RealVec )
    return dot(p,p)
end

"""
    Enf( m::Real, ps::Real )

Kinetic energy function. Returns ``E=\\sqrt{m^2+p^2}``, given ``m``
and ``\\vec{p}``.
"""
function Enf( m::Real, ps::Real )
    return sqrt(m^2+ps)
end

"""
    rf( qa::RealVec , qb::RealVec )

Euclidean distance function. Returns ``r_{ab}=\\sqrt{(qa-qb)⋅(qa-qb)}``,
given particle positions ``\\vec{q}_a`` and ``\\vec{q}_b``.
"""
function rf( qa::RealVec , qb::RealVec )
    return norm(qa-qb)
end

"""
    nabf( qa::RealVec , qb::RealVec )

Separation unit vector. Returns 
``\\vec{n}_{ab}=(\\vec{q}_a-\\vec{q}_b)/\\sqrt{r_{ab}}``,
given particle positions ``\\vec{q}_a`` and ``\\vec{q}_b``.
"""
function nabf( qa::RealVec , qb::RealVec )
    return (qa-qb)/rf(qa,qb)
end

"""
    ybaf( mb::Real, qb::RealVec , qa::RealVec , pb::RealVec , pa::RealVec )

Returns ``y_{ba}=\\sqrt{m_b^2+(\\vec{n}_{ab}⋅\\vec{p}_b)^2}/E_b``, where ``E_b`` is
Kinetic energy for particle ``b``. Inputs are particle ``b`` mass ``m_b``,
particle positions ``\\vec{q}_b`` and ``\\vec{q}_a``, and particle momenta
``\\vec{p}_b`` and ``\\vec{p}_a``.
"""
function ybaf( mb::Real, qb::RealVec , qa::RealVec , pb::RealVec , pa::RealVec )
    return sqrt(mb^2+dot(nabf(qb,qa),pb)^2)/Enf(mb,psf(pb))
end

"""
    Θabf( qa::RealVec , qb::RealVec , pa::RealVec )

Returns ``Θ_{ba}=\\vec{p}_a⋅\\vec{n}_{ab}``, where ``\\vec{p}_b`` is
momentum for particle ``b``, and ``\\vec{n}_{ab}`` is unit separation 
vector.
"""
function Θabf( qa::RealVec , qb::RealVec , pa::RealVec )
    return dot(pa,nabf(qb,qa))
end

"""
    Ξabf( pa::RealVec, pb::RealVec )

Returns ``Ξ_{ba}=\\vec{p}_a⋅\\vec{p}_{b}``, where ``\\vec{p}_b`` and 
``\\vec{p}_a`` are momenta for particles ``a`` and ``b``.
"""
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
                    H3 += ( 2*(-2*(Ξab*Θba)^2 - 2*Θab*Θba*Ξab*psb - (Θab*psb)^2 + psb*Ξab^2)/Enb^2 + 
                            2*(psa*Θba^2 - (Θab*Θba)^2 + 2*Θab*Θba*Ξab - Ξab^2 + psb*Θab^2) +
                            yba*(3*psa*Θba^2 - (Θab*Θba)^2 + 8*Θab*Θba*Ξab - psa*psb + 3*psb*Θab^2)
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

    return H0+H1+H2+H3

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

# Gradient of Milky Way potentials
function dMilkyWay(d::Int, m::RealVec, Z::RealVec, origin_x=-1.708859462494220E+17, origin_y=0, origin_z =4.346342845091530E+14, r_sun=8.4, Mb=409, Md=2856, Mh=1018, b_b=0.23, a_d=4.22, b_d=0.292, a_h=2.562)
    # (origin_x, origin_y, origin_z) is the location (in units of M) in the galactocentric frame where the simulation's origin is placed
    # it defaults to the location of the Sun as given by astropy and 2019 data from https://arxiv.org/abs/1904.05721
    
    # Milky Way Model I (Irrgang et al):
    #   r_sun = 8.4 kpc       radius of sun's orbit -- not used
    #   Mb = 409 M_gal        mass bulge
    #   Md = 2856 M_gal       mass disk
    #   Mh = 1018 M_gal       mass halo
    #   b_b = 0.23 kpc        length scale bulge
    #   a_d = 4.22 kpc        length scale disk #1
    #   b_d = 0.292 kpc       length scale disk #2
    #   a_h = 2.562 kpc       length scale halo
    #   Λ = 200 kpc           halo cutoff parameter
    #   γ = 2                 free parameter
    
    tpfl = typeof(Z[1])
    n = length(m)

    # the Milky Way potentials are in 3-dimensionsal Cartesian coordinates
    # ensure d is 3, otherwise return all zeros (effectively ignoring Milky Way potentials)
    if d != 3 
        return zeros(tpfl, 6*n)
    end

    dΦB = zeros(tpfl, 6*n)
    dΦD = zeros(tpfl, 6*n)
    dΦH = zeros(tpfl, 6*n)

    # Z consists of q's for all particles first, then p's for all particles
    # dΦ vectors are gradients wrt Z, therefore derivs wrt all of the x, y, z's are first, then wrt all of the px, py, pz's
    # Since MW Φ's contain no p's, the second half of each gradient vector will be zero (whole vector already initialized to zeros)
    for i in 1:n
        # index for start of q's of particle i in Z
        part_index = 1 + 3*(i-1)  # since d == 3
        # indices for x, y, z positions of particle i within Z, which also are indices within dΦ for x, y, z derivs for particle i
        xind = part_index
        yind = part_index + 1
        zind = part_index + 2
        pos = [ Z[xind] + origin_x, Z[yind] + origin_y, Z[zind] + origin_z]
        
        # bulge
        bulge_denom = ( dot(pos,pos) + b_b^2)^1.5
        dΦB[xind] += Mb * pos[1] / bulge_denom
        dΦB[yind] += Mb * pos[2] / bulge_denom
        dΦB[zind] += Mb * pos[3] / bulge_denom

        # disc
        rsq = pos[1]^2 + pos[2]^2
        disc_denom = ( rsq + (a_d + sqrt(pos[3]^2 + b_d^2))^2)^1.5
        dΦD[xind] += Md * pos[1] / disc_denom
        dΦD[yind] += Md * pos[2] / disc_denom
        dΦD[zind] += Md * pos[3] * (a_d + sqrt(pos[3]^2 + b_d^2)) / ( disc_denom * sqrt(pos[3]^2 + b_d^2))

        # halo
        R = sqrt(dot(pos,pos))
        halo_coeff = Mh / a_h^2 * (1 + R/a_h)^(-1) / R
        dΦH[xind] += halo_coeff * pos[1]
        dΦH[yind] += halo_coeff * pos[2]
        dΦH[zind] += halo_coeff * pos[3]
    end

    return dΦB + dΦD + dΦH
end

# Gradient of the Hamiltonian function
function dH( d::Int , m::RealVec , Z::RealVec )
    return ForwardDiff.gradient(x -> H(d, m, x), Z) #+ dH3m0( d , m , Z )
end

# Gradient of the Hamiltonian function plus gradient of Milky Way potential
function dH_plus_MW(d::Int, m::RealVec, Z::RealVec, origin_x=-1.708859462494220E+17, origin_y=0, origin_z =4.346342845091530E+14)
    # (origin_x, origin_y, origin_z) is the location (in units of M) in the galactocentric frame where the simulation's origin is placed
    # it defaults to the location of the Sun as given by astropy and 2019 data from https://arxiv.org/abs/1904.05721

    return ForwardDiff.gradient(x -> H(d, m, x), Z) + dMilkyWay(d, m, Z, origin_x, origin_y, origin_z, 1.7552537847E+17, 9.5091683066E+09, 6.6401429544E+10, 2.3668296665E+10, 4.8060520294E+15, 8.8180606801E+16, 6.1015964896E+15, 5.3535240432E+16)
end

## Symplectic operator: Maps output of dH to time derivative of phase space variables
#function Jsympl( Zarg::RealVec )
#    tpfl=typeof(Zarg[1])
#    n2 = length(Zarg)
#    Z = zeros(tpfl,n2)

#    if iseven(n2)
#        n = Int(round(n2 / 2, digits=0))
#        for i=1:n
#            Z[i]    = Zarg[n+i] 
#            Z[n+i]  = - Zarg[i] 
#        end
#
#        return Z
#    else
#        return Z
#    end
#end

# Right hand side of Hamilton's equations
FHE = (d,m,z)->Jsympl(dH(d,m,z))

#-----------------------------------------------------------------------
#   SIMPLE ADAPTIVE TIMESTEPPING
#-----------------------------------------------------------------------
"""
    tcour( dt::Real , Z::RealVec , Zdot::RealVec , C = 0.001 , d=3 )

This function implements a simple adaptive timestepping function 
inspired by the Courant-Friedrichs Condition. 

Returns a new timestep that satisfies the CFL condition.
If the current timestep already satisfies CFL condition, the timestep is returned unchanged.

# Arguments
- `dt::Real` is the current timestep
- `Z::RealVec` is the current state of the system
- `Zdot::RealVec` is the time derivative of Z
- `C` is the Courant factor
- `d` is the number of dimensions

"""
function tcour( dt::Real , Z::RealVec , Zdot::RealVec , C = 0.001 , d=3 )
    tpfl = typeof(Z[1])
    n2 = length(Z)
    rs = tpfl(1)
    vs = tpfl(1)

    ndof = Int(round(n2 / 2, digits=0))
    n    = Int(round(ndof/d, digits=0))

    if iseven(n2) && n*d==ndof
        # calculate separation and relative velocity between first two particles
        # and store these as starting minimum values
        qa  = Z2q(n,d,1,Z)
        qb  = Z2q(n,d,2,Z)
        va  = Z2q(n,d,1,Zdot)
        vb  = Z2q(n,d,2,Zdot)
        vab = norm(va-vb)
        rab = rf(qa,qb)
        rs  = rab  # smallest separation among particle pairs
        vs  = vab       # relative velocity for particle pair with smallest separation
        for a=1:n
            qa = Z2q(n,d,a,Z)    
            va = Z2q(n,d,a,Zdot)
            for b=2:n
                if b!=a
                    qb = Z2q(n,d,b,Z)
                    vb = Z2q(n,d,b,Zdot)
                    rab = rf(qa,qb)
                    vab = norm(va-vb)
                    if rab<=rs  
                        vs = vab
                        rs = rab
                    end
                end
            end
        end
        Δt = C * rs/vs
        # return \Delta t if it is smaller than the maximum timestep, otherwise return maximum timestep

        # println("time step = " * string(Δt))
        return Δt
    else
        error("in tcour delta t is not changing bc ( iseven(n2) && n*d==ndof ) returned false")
        return dt
    end
end  # End tcour
