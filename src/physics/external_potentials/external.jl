#-----------------------------------------------------------------------
#   MILKY WAY POTENTIAL GRADIENTS
#-----------------------------------------------------------------------

# Gradient of Milky Way potentials
function ΦMilkyWay(Z::RealVec, origin_x=Double64(-1.708859462494220E+17), origin_y=Double64(0), origin_z =Double64(4.346342845091530E+14), r_sun=Double64(8.4), Mb=Double64(409), Md=Double64(2856), Mh=Double64(1018), b_b=Double64(0.23), a_d=Double64(4.22), b_d=Double64(0.292), a_h=Double64(2.562), Λ=Double64(100), γ=Double64(2.02),d::Int=3)
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
    n = round(Int,Float64(length(Z)/(2*d)))

    PHI = zero(tpfl)

    if d == 3 
        for i in 1:n
            qi = Z2q( n , d , i , Z )
    
            Q  = qi .+ [origin_x, origin_y, origin_z]
    
            R     = norm(Q)
            r     = sqrt(Q[1]^2 + Q[2]^2)
            z     = Q[3]
    
            PHI_b = -Mb/sqrt(R^2+b_b^2)
    
            PHI_d = -Md/sqrt(r^2 + (a_d + sqrt(z^2 + b_d^2))^2)
    
            PHI_h = (Mh/ah) * 
                    ( 
                     (1/(γ-1))*log((1 + (R/ah)^(γ-1))/(1 + (Λ/ah)^(γ-1))) - 
                     (Λ/ah)^(γ-1) / (1 + (Λ/ah)^(γ-1)) 
                    )
    
            PHI += PHI_b + PHI_d + PHI_h
        end
        return PHI
    else 
        return zero(tpfl)
    end
end #-------------------------------------------------------------------

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
