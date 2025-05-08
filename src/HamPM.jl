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

Kinetic energy function. Returns ``E=\\sqrt{m^2+p^2}`` for massive particles (m>0),
and ``E=\\sqrt{p^2}=|p|`` for massless particles (m=0), given ``m`` and ``\\vec{p}``.
"""
function Enf( m::Real, ps::Real )
    if m > zero(m)
        return sqrt(m^2+ps)
    else
        return sqrt(ps)
    end
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
    return sqrt(mb^2+Θabf(qb,qa,pb)^2)/Enf(mb,psf(pb))
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

                if m[b] != zero(tpfl) && yba != zero(tpfl) && yba != one(tpfl)
                    H3 += ( 2*(-2*(Ξab*Θba)^2 - 2*Θab*Θba*Ξab*psb - (Θab*psb)^2 + psb*Ξab^2)/Enb^2 + 
                            2*(psa*Θba^2 - (Θab*Θba)^2 + 2*Θab*Θba*Ξab - Ξab^2 + psb*Θab^2) +
                            yba*(3*psa*Θba^2 - (Θab*Θba)^2 + 8*Θab*Θba*Ξab - psa*psb + 3*psb*Θab^2)
                            ) / (4*Ena*Enb*rab*yba*(yba+1)^2)
                else
                    H3 += (3*Enb*psb*Θab^2 - Enb*psa*(psb - 3*Θba^2) + Enb*Θab*Θba*(-(Θab*Θba) + 8*Ξab) + 2*psa*psb*abs(Θba) - 2*psb*Θab^2*abs(Θba) - 4*Ξab^2*abs(Θba))/(4*Ena*rab*(Enb + abs(Θba))^2)
                end
            end
        end
    end

    return H0+H1+H2+H3

end

#-----------------------------------------------------------------------
#   POST-MINKOWSKIAN HAMILTONIAN - MASSLESS PARTICLE IMPLEMENTATION
#-----------------------------------------------------------------------

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

"""
    compute_dy_massless(Ena, dTha, dEa, Tha)

Compute the y-parameter derivative for massless particles.
This matches the C implementation's calculation.
"""
function compute_dy_massless(Ena, dTha, dEa, Tha)
    return sign(Tha)*(Ena*dTha - dEa*Tha)/(Ena*Ena)
end

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
                # First loop over a (c fixed)
                for a=1:n
                    if a != c
                        qa = Z2q(n,d,a,Z)
                        qc = Z2q(n,d,c,Z)
                        pa = Z2p(n,d,a,Z)
                        pc = Z2p(n,d,c,Z)

                        rac = sqrt(sum((qa .- qc).^2))
                        psa = sqrt(sum(pa.^2))
                        psc = sqrt(sum(pc.^2))
                        Ena = sqrt(m[a]^2 + psa^2)
                        Enc = psc  # For massless particle c

                        # Calculate Theta and Xi
                        Θac = sum((qa .- qc).*pc)/rac
                        Θca = -sum((qa .- qc).*pa)/rac
                        Ξac = sum(pa.*pc)

                        # Compute derivatives
                        dr_ac = (qc[i]-qa[i])/rac
                        dTh_ac = (pa[i]-Θac*dr_ac)/rac
                        dXi_ac = zero(tpfl)

                        # Use proper y-parameter derivative for massless case
                        dy_ac = compute_dy_massless(Ena, dTh_ac, zero(tpfl), Θac)

                        # Add contribution from a terms
                        dH3[i+d*(c-1)] += compute_dH3_term(Ena, Enc, dy_ac, rac, sign(Θca), psc, Θac, Θca, Ξac, psa, Ena, dr_ac, dTh_ac, dXi_ac)
                    end
                end

                # Second loop over b (c fixed)
                for b=1:n
                    if b != c
                        qb = Z2q(n,d,b,Z)
                        qc = Z2q(n,d,c,Z)
                        pb = Z2p(n,d,b,Z)
                        pc = Z2p(n,d,c,Z)

                        rbc = sqrt(sum((qb .- qc).^2))
                        psb = sqrt(sum(pb.^2))
                        psc = sqrt(sum(pc.^2))
                        Enb = sqrt(m[b]^2 + psb^2)
                        Enc = psc  # For massless particle c

                        # Calculate Theta and Xi
                        Θbc = sum((qb .- qc).*pc)/rbc
                        Θcb = -sum((qb .- qc).*pb)/rbc
                        Ξbc = sum(pb.*pc)

                        # Compute derivatives
                        dr_bc = (qc[i]-qb[i])/rbc
                        dTh_bc = (pb[i]-Θbc*dr_bc)/rbc
                        dXi_bc = zero(tpfl)

                        # Use proper y-parameter derivative for massless case
                        dy_bc = compute_dy_massless(Enb, dTh_bc, zero(tpfl), Θbc)

                        # Add contribution from b terms
                        dH3[i+d*(c-1)] += compute_dH3_term(Enb, Enc, dy_bc, rbc, sign(Θcb), psc, Θbc, Θcb, Ξbc, psb, Enb, dr_bc, dTh_bc, dXi_bc)
                    end
                end
            end

            # Momentum derivatives (zindex >= 3)
            for i=1:d
                # Calculate momentum derivatives once per dimension
                dps_c = 2*pc[i]
                dE_c = dps_c/(2*Enc)

                # First loop over a (c fixed)
                for a=1:n
                    if a != c
                        qa = Z2q(n,d,a,Z)
                        qc = Z2q(n,d,c,Z)
                        pa = Z2p(n,d,a,Z)
                        pc = Z2p(n,d,c,Z)

                        rac = sqrt(sum((qa .- qc).^2))
                        psa = sqrt(sum(pa.^2))
                        psc = sqrt(sum(pc.^2))
                        Ena = sqrt(m[a]^2 + psa^2)
                        Enc = psc  # For massless particle c

                        # Calculate Theta and Xi
                        Θac = sum((qa .- qc).*pc)/rac
                        Θca = -sum((qa .- qc).*pa)/rac
                        Ξac = sum(pa.*pc)

                        # Compute derivatives
                        dTh_ac = (pa[i]-Θac*zero(tpfl))/rac
                        dXi_ac = zero(tpfl)

                        # Use proper y-parameter derivative for massless case
                        dy_ac = compute_dy_massless(Ena, dTh_ac, zero(tpfl), Θac)

                        # Add contribution from a terms
                        dH3[i+d*(n+c-1)] += compute_dH3_term(Ena, Enc, dy_ac, rac, sign(Θca), psc, Θac, Θca, Ξac, psa, dE_c, zero(tpfl), dTh_ac, dXi_ac)
                    end
                end

                # Second loop over b (c fixed)
                for b=1:n
                    if b != c
                        qb = Z2q(n,d,b,Z)
                        qc = Z2q(n,d,c,Z)
                        pb = Z2p(n,d,b,Z)
                        pc = Z2p(n,d,c,Z)

                        rbc = sqrt(sum((qb .- qc).^2))
                        psb = sqrt(sum(pb.^2))
                        psc = sqrt(sum(pc.^2))
                        Enb = sqrt(m[b]^2 + psb^2)
                        Enc = psc  # For massless particle c

                        # Calculate Theta and Xi
                        Θbc = sum((qb .- qc).*pc)/rbc
                        Θcb = -sum((qb .- qc).*pb)/rbc
                        Ξbc = sum(pb.*pc)

                        # Compute derivatives
                        dTh_bc = (pb[i]-Θbc*zero(tpfl))/rbc
                        dXi_bc = zero(tpfl)

                        # Use proper y-parameter derivative for massless case
                        dy_bc = compute_dy_massless(Enb, dTh_bc, zero(tpfl), Θbc)

                        # Add contribution from b terms
                        dH3[i+d*(n+c-1)] += compute_dH3_term(Enb, Enc, dy_bc, rbc, sign(Θcb), psc, Θbc, Θcb, Ξbc, psb, dE_c, zero(tpfl), dTh_bc, dXi_bc)
                    end
                end
            end
        end
    end

    return dH3
end

#-----------------------------------------------------------------------
#   MIKLY WAY POTENTIAL GRADIENTS
#-----------------------------------------------------------------------

# Gradient of Milky Way potentials
function dMilkyWay(d::Int, m::RealVec, Z::RealVec, origin_x=Double64(-1.708859462494220E+17), origin_y=Double64(0), origin_z =Double64(4.346342845091530E+14), r_sun=Double64(8.4), Mb=Double64(409), Md=Double64(2856), Mh=Double64(1018), b_b=Double64(0.23), a_d=Double64(4.22), b_d=Double64(0.292), a_h=Double64(2.562))
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
        # index in Z for start of q's of particle i
        part_index = 1 + 3*(i-1)  # since d == 3
        # indices for x, y, z positions of particle i within Z, which also are indices within dΦ for x, y, z derivs for particle i
        xind = part_index
        yind = part_index + 1
        zind = part_index + 2

        pos = [ Z[xind] + origin_x, Z[yind] + origin_y, Z[zind] + origin_z]
        
        # bulge
        bulge_denom = ( dot(pos,pos) + b_b^2)^(Double64(1.5))
        dΦB[xind] += m[i] * Mb * pos[1] / bulge_denom
        dΦB[yind] += m[i] * Mb * pos[2] / bulge_denom
        dΦB[zind] += m[i] * Mb * pos[3] / bulge_denom

        # disc
        rsq = pos[1]^2 + pos[2]^2
        disc_denom = ( rsq + (a_d + sqrt(pos[3]^2 + b_d^2))^2)^(Double64(1.5))
        dΦD[xind] += m[i] * Md * pos[1] / disc_denom
        dΦD[yind] += m[i] * Md * pos[2] / disc_denom
        dΦD[zind] += m[i] * Md * pos[3] * (a_d + sqrt(pos[3]^2 + b_d^2)) / (disc_denom * sqrt(pos[3]^2 + b_d^2))

        # halo
        R = sqrt(dot(pos,pos))
        halo_coeff = m[i] * Mh / a_h^2 * (1 + R / a_h)^(-1) / R
        dΦH[xind] += halo_coeff * pos[1]
        dΦH[yind] += halo_coeff * pos[2]
        dΦH[zind] += halo_coeff * pos[3]
    end

    dΦ = dΦB + dΦD + dΦH
    return dΦ
end

# Computes gradient of Sun's potential at position of each particle.  Assumes Sun at origin with M = 1
function dSun(d::Int, m::RealVec, Z::RealVec)

    tpfl = typeof(Z[1])
    n = length(m)

    # ensure d is 3, otherwise return all zeros (effectively ignoring potentials)
    if d != 3
        return zeros(tpfl, 6 * n)
    end

    dS = zeros(tpfl, 6 * n)

    for i in 1:n
        # index in Z for start of q's of particle i
        part_index = 1 + 3 * (i - 1)  # since d == 3
        # indices for x, y, z positions of particle i within Z, which also are indices within dS for x, y, z derivs for particle i
        xind = part_index
        yind = part_index + 1
        zind = part_index + 2

        pos = [Z[xind], Z[yind], Z[zind]]
        r_cubed = norm(pos) * norm(pos) * norm(pos)
        
        # G == M == 1
        dS[xind] = m[i] * pos[1] / r_cubed
        dS[yind] = m[i] * pos[2] / r_cubed
        dS[zind] = m[i] * pos[3] / r_cubed
    end

    return dS
end

# Gradient of the Hamiltonian function
function dH( d::Int , m::RealVec , Z::RealVec )
    return ForwardDiff.gradient(x -> H(d, m, x), Z) + dH3m0( d , m , Z )
end

# Gradient of the Hamiltonian function plus gradient of Milky Way potential
function dH_plus_MW(d::Int, m::RealVec, Z::RealVec, origin_x=Double64(-1.708859462494220E+17), origin_y=Double64(0), origin_z=Double64(4.346342845091530E+14))
    # (origin_x, origin_y, origin_z) is the location (in units of M) in the galactocentric frame where the simulation's origin is placed
    # it defaults to the location of the Sun as given by astropy and 2019 data from https://arxiv.org/abs/1904.05721

    grad = ForwardDiff.gradient(x -> H(d, m, x), Z) 
    dMW = dMilkyWay(d, m, Z, origin_x, origin_y, origin_z, Double64(1.7552537847E+17), Double64(9.5091683066E+09), Double64(6.6401429544E+10), Double64(2.3668296665E+10), Double64(4.8060520294E+15), Double64(8.8180606801E+16), Double64(6.1015964896E+15), Double64(5.3535240432E+16))
    
    return grad + dMW    
end

function dH_plus_Sun(d::Int, m::RealVec, Z::RealVec)

    grad = ForwardDiff.gradient(x -> H(d, m, x), Z) 
    dS = dSun(d, m, Z)
    
    return grad + dS
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
# FHE = (d,m,z)->Jsympl(dH(d,m,z))

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

function dt_lin_int(dt::Real , Z::RealVec , Zdot::RealVec , start_dist::Real , min_dist::Real , max_dt::Real , min_dt::Real, d = 3)
"""
    This function implements a simple adaptive timestepping function by linear interpolation.
    When distance between first two bodies is equal to or greater than start_dist, it uses the max_dt
    When distance is equal to or less than min_dist, it uses the min_dt
    Between the two, it interpolates linearly
"""
    # find distance between first two bodies

    tpfl = typeof(Z[1])
    n2 = length(Z)

    ndof = Int(round(n2 / 2, digits=0))
    n    = Int(round(ndof/d, digits=0))
    if iseven(n2) && n*d==ndof
        qa  = Z2q(n,d,1,Z)
        qb  = Z2q(n,d,2,Z)
        rab = rf(qa,qb)
    else
        error("in dt_lin_int delta t is not changing bc ( iseven(n2) && n*d==ndof ) returned false")
        return dt
    end

    if rab >= start_dist
        dt = max_dt
    elseif rab <= min_dist
        dt = min_dt
    else
        # interpolate linearly
        dt = abs(start_dist - rab) / abs(start_dist - min_dist) * min_dt + abs(rab - min_dist) / abs(start_dist - min_dist) * max_dt
    end

    # double-check that we haven't gone outside range of [min_dt, max_dt]
    if dt > max_dt
        dt = max_dt
    end
    if dt <= min_dt
        dt = min_dt
    end

    return dt

end
