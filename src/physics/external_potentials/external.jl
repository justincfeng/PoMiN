#-----------------------------------------------------------------------
#   MILKY WAY POTENTIAL GRADIENTS
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
