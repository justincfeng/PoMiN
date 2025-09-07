#-----------------------------------------------------------------------
#   POTENTIAL ENERGY CONSTRUCTOR FUNCTION
#-----------------------------------------------------------------------

"""
    UConstructor(PHI::Function,m::RealVec,xo::RealVec,
                              sc::Real=1.0,d::Int=3)

This function takes a potential function 'PHI(x)' and returns a new
potential energy function that is the sum of 'PHI(x)' evaluated at each
particle position 'x[i]' multiplied by the mass of the particle 'm[i]'.
It is used to construct external potentials.
"""
function UConstructor(PHI::Function,m::RealVec,xo::RealVec,
                              sc::Real=1.0,d::Int=3)
    tpfl = typeof(m[1])
    n = length(m)

    function U(Z::RealVec)
        V = zero(eltype(Z))  # Use eltype(Z) for ForwardDiff compatibility
        for i in 1:n
            x = (sc .* Z2q( n , d , i , Z )) + xo  # Remove type conversion for ForwardDiff compatibility
            V += m[i]*PHI(x)
        end
        return V
    end
    return U
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   MILKY WAY POTENTIAL
#-----------------------------------------------------------------------

# Milky Way potential for a single particle
"""
    ΦMilkyWay( tpfl::Type=Double64 ,
                origin_x=tpfl(-1.708859462494220E+17), origin_y=tpfl(0), 
                origin_z =tpfl(4.346342845091530E+14), r_sun=tpfl(8.4), 
                Mb=tpfl(409), Md=tpfl(2856), Mh=tpfl(1018), 
                b_b=tpfl(0.23), a_d=tpfl(4.22), b_d=tpfl(0.292), 
                a_h=tpfl(2.562), Λ=tpfl(100), γ=tpfl(2.02),d::Int=3 )

This function returns a potential function for a single particle under
the influence of the gravitational potential of the Milky Way. The 
assumed length units for the argument are in terms of solar masses 
in geometric units (G=c=1).
"""
function ΦMilkyWay( tpfl::Type=Double64 ,
                origin_x=tpfl(0.0), origin_y=tpfl(0.0), 
                origin_z=tpfl(0.0), r_sun=tpfl(8.4), 
                Mb=tpfl(409), Md=tpfl(2856), Mh=tpfl(1018), 
                b_b=tpfl(0.23), a_d=tpfl(4.22), b_d=tpfl(0.292), 
                a_h=tpfl(2.562), Λ=tpfl(200), γ=tpfl(2.0),d::Int=3 )
    # Corrected Milky Way Model I implementation based on Irrgang et al. (arXiv:1211.4353)
    # Default origin at (0,0,0) for galactocentric coordinates
    # Original solar offset: origin_x=-1.708859462494220E+17, origin_z=4.346342845091530E+14
    
    # Milky Way Model I parameters (Irrgang et al., Table 1):
    #   r_sun = 8.4 kpc       radius of sun's orbit
    #   Mb = 409 M_gal        bulge mass parameter  
    #   Md = 2856 M_gal       disk mass parameter
    #   Mh = 1018 M_gal       halo mass parameter
    #   b_b = 0.23 kpc        bulge scale length
    #   a_d = 4.22 kpc        disk scale length #1
    #   b_d = 0.292 kpc       disk scale length #2
    #   a_h = 2.562 kpc       halo scale length
    #   Λ = 200 kpc           halo cutoff parameter
    #   γ = 2                 halo exponent (fixed)

    # Unit conversions between paper and PoMiN systems
    Mgal = tpfl(2.325e7)                    # M_gal in M_sun
    Msol2m = tpfl(1476.67)                  # 1 M_sun in meters  
    kpc_m = tpfl(3.0857e19)                 # 1 kpc in meters
    kpc = kpc_m / Msol2m                    # 1 kpc in solar mass units
    
    # Convert potential units: 100 km²/s² to PoMiN units (c=1)
    c_mks = tpfl(299792458.0)
    potential_unit_conversion = tpfl(1e8) / (c_mks * Msol2m)^2

    if d == 3
        ν1 = tpfl(1)
        function Φ(x::RealVec)
            # Apply origin offset
            X = x .+ [origin_x, origin_y, origin_z]
            
            # Convert to kpc for potential calculation
            x_kpc = X[1] / kpc
            y_kpc = X[2] / kpc  
            z_kpc = X[3] / kpc
            
            # Coordinate calculations
            R = sqrt(x_kpc^2 + y_kpc^2 + z_kpc^2)  # Spherical radius
            r = sqrt(x_kpc^2 + y_kpc^2)            # Cylindrical radius
            z = z_kpc                               # Height
            
            # Bulge potential (Miyamoto-Nagai spherical)
            PHI_b = -Mb / sqrt(R^2 + b_b^2)
            
            # Disk potential (Miyamoto-Nagai disk)
            PHI_d = -Md / sqrt(r^2 + (a_d + sqrt(z^2 + b_d^2))^2)
            
            # Halo potential (Allen & Santillan with γ=2)
            if R < Λ
                term1 = (ν1/(γ-ν1)) * log((ν1 + (R/a_h)^(γ-ν1)) / (ν1 + (Λ/a_h)^(γ-ν1)))
                term2 = (Λ/a_h)^(γ-ν1) / (ν1 + (Λ/a_h)^(γ-ν1))
                PHI_h = (Mh/a_h) * (term1 - term2)
            else
                # Outside cutoff radius
                PHI_h = -(Mh/R) * ((Λ/a_h)^γ / (1 + (Λ/a_h)^(γ-1)))
            end
            
            # Total potential in paper units (100 km²/s²)
            PHI_total_paper = PHI_b + PHI_d + PHI_h
            
            # Convert to PoMiN units
            PHI_total = PHI_total_paper * potential_unit_conversion
            
            return PHI_total
        end
        return Φ
    else 
        return x->zero(eltype(x))
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   CENTRAL MASS POTENTIAL
#-----------------------------------------------------------------------
"""
    ΦCentralMass( xo::RealVec , tpfl::Type=Double64 , 
                  M::Real=1.0 , G::Real=1.0 , d::Int=3 )

This function returns a potential function for a single particle under
the influence of the gravitational potential of a central mass.
"""
function ΦCentralMass(xo::RealVec , tpfl::Type=Double64 , 
                      M::Real=1.0 , G::Real=1.0 , d::Int=3 )
    return x->-G*M/norm(x-xo)  # Remove type conversion to allow ForwardDiff compatibility
end #-------------------------------------------------------------------

