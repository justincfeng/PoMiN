#-----------------------------------------------------------------------
#   POST-MINKOWSKIAN HAMILTONIAN - MASSLESS PARTICLE IMPLEMENTATION
#-----------------------------------------------------------------------

"""
    dH3m0( d::Int , m::RealVec , Z::RealVec )

Derivative of H3 for massless particles. This function handles the special case
where particle c has zero mass.
"""
function dH3m0( d::Int , m::RealVec , Z::RealVec )
    tpfl = typeof(Z[1])
    n = length(m)
    nn = 2*n*d

    o = zero(tpfl)

    # Initialize arrays for derivatives
    dps = zeros(tpfl,n)
    dE = zeros(tpfl,n)
    dr = [zeros(tpfl,n) for _ = 1:n]
    dy = [zeros(tpfl,n) for _ = 1:n]
    dTh = [zeros(tpfl,n) for _ = 1:n]
    dXi = [zeros(tpfl,n) for _ = 1:n]

    dH3 = zeros(tpfl,nn)

    # Loop over particles
    for c=1:n
        # Only proceed if particle c is massless
        if m[c] == o
            # Position derivatives (zindex < 3)
            for i=1:d
                # Loop over other particles for cross-terms
                for a=1:n
                    if a != c
                        qa = Z2q(n,d,a,Z)
                        qc = Z2q(n,d,c,Z)
                        pa = Z2p(n,d,a,Z)
                        pc = Z2p(n,d,c,Z)
                        
                        # Calculate basic quantities
                        psa = psf(pa)
                        psc = psf(pc)
                        Ena = Enf(m[a],psa)
                        Enc = Enf(m[c],psc)
                        rac = rf(qa,qc)
                        
                        # Calculate Theta and Xi
                        Θac = Θabf(qa,qc,pa)
                        Θca = Θabf(qc,qa,pc)
                        Ξac = Ξabf(pa,pc)
                        
                        # Calculate y for massless case
                        yca = sign(Θca)
                        
                        # Position derivatives
                        dr[a][c] = (qc[i]-qa[i])/rac
                        dTh[a][c] = (pa[i]-Θac*dr[a][c])/rac
                        dXi[a][c] = o
                        
                        # Special handling for massless y derivative
                        dyac = sign(Θac)*(Ena*dTh[a][c] - dE[a]*Θac)/(Ena^2)
                        
                        # Momentum derivatives
                        if i == 1  # Only calculate once per dimension
                            dps[c] = 2*pc[i]
                            dE[c] = dps[c]/(2*Enc)
                        end
                        
                        # Update dH3 for both position and momentum components
                        dH3[d*(c-1)+i] -= compute_dH3_term(Ena, Enc, dyac, rac, yca, psc, Θac, Θca, Ξac, psa, dE, dr, dTh, dXi)
                        dH3[d*(c-1+n)+i] -= compute_dH3_term(Ena, Enc, dyac, rac, yca, psc, Θac, Θca, Ξac, psa, dE, dr, dTh, dXi)
                    end
                end
            end
        end
    end

    return dH3
end

"""
    compute_dH3_term(Ena, Enc, dyac, rac, yca, psc, Θac, Θca, Ξac, psa, dE, dr, dTh, dXi)

Helper function to compute the H3 derivative terms for massless particles.
This matches the C implementation's calculation for the massless case.
"""
function compute_dH3_term(Ena, Enc, dyac, rac, yca, psc, Θac, Θca, Ξac, psa, dE, dr, dTh, dXi)
    return (1/4)*(-((2*Ena*Enc*dyac*rac*yca*(2*psc*psc*Θac*Θac + 4*psc*Θac*Θca*Ξac - 
           2*(psc - 2*Θca*Θca)*Ξac*Ξac + Enc*Enc*(2*(-(Θac*Θca) + Ξac)*(-(Θac*Θca) + Ξac) + 
           Θac*Θca*(Θac*Θca - 8*Ξac)*yca - psa*Θca*Θca*(2 + 3*yca) + 
           psc*(psa*yca - Θac*Θac*(2 + 3*yca))))))/(Ena*Ena*Enc*Enc*Enc*Enc*rac*rac*yca*yca*(1 + yca)*(1 + yca)*(1 + yca)))
end
