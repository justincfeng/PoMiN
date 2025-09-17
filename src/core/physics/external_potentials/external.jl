#-----------------------------------------------------------------------
#   MILKY WAY POTENTIAL
#-----------------------------------------------------------------------

# Milky Way potential for a single particle
"""

    ΦMilkyWay( tpfl::Type=Double64 , x0::RealVec=zeros(tpfl,3)     ,
                p::Tuple = ( tpfl(409), tpfl(2856), tpfl(1018)   , 
                             tpfl(0.23), tpfl(4.22), tpfl(0.292) , 
                             tpfl(2.562), tpfl(200), tpfl(2) )     )

This function returns a potential function for a single particle under
the influence of the gravitational potential of the Milky Way. The 
assumed length units for the argument are in terms of solar masses 
in geometric units (G=c=1). This function uses the Milky Way Model I 
from Irrgang et al. (arXiv:1211.4353). The arguments are as follows:

    tpfl::Type  : The floating point type of the potential function
    x0::RealVec : The origin of the potential function
    p::Tuple    : The parameters of the potential function

The parameters `p` of the potential function are as follows:

    p[1] : Mb    : The bulge mass parameter
    p[2] : Md    : The disk mass parameter
    p[3] : Mh    : The halo mass parameter
    p[4] : b_b   : The bulge scale length
    p[5] : a_d   : The disk scale length #1
    p[6] : b_d   : The disk scale length #2
    p[7] : a_h   : The halo scale length
    p[8] : Λ     : The halo cutoff parameter
    p[9] : γ     : The halo exponent (fixed)

The assumed units for p are in kpc and in galactic mass Mgal, which are
then converted to geometric units (G=c=1).

"""
function ΦMilkyWay( tpfl::Type=Double64 , x0::RealVec=zeros(tpfl,3)    ,
                    p::Tuple = ( tpfl(409), tpfl(2856), tpfl(1018)   , 
                                 tpfl(0.23), tpfl(4.22), tpfl(0.292) , 
                                 tpfl(2.562), tpfl(200), tpfl(2) )     )
    # FOR REFERENCE:
    # solar offset: [-1.708859462494220E+17, 0.0, 4.346342845091530E+14]

    Mb_mg, Md_mg, Mh_mg, b_b_kpc, a_d_kpc, b_d_kpc, a_h_kpc, Λ_kpc, γ = p

    # Unit conversions between paper and PoMiN systems
    Mgal = tpfl(2.325e7)                    # M_gal in M_sun
    Msol2m = tpfl(1476.67)                  # 1 M_sun in meters  
    kpc_m = tpfl(3.0857e19)                 # 1 kpc in meters
    kpc = kpc_m / Msol2m                    # 1 kpc in solar mass units

    Mb , Md , Mh = Mb_mg * Mgal , Md_mg * Mgal , Mh_mg * Mgal
    a_d , a_h = a_d_kpc * kpc , a_h_kpc * kpc
    b_b , b_d = b_b_kpc * kpc , b_d_kpc * kpc
    Λ = Λ_kpc * kpc
    ν1 = tpfl(1)

    function Φ(x::RealVec)
        # Apply origin offset
        X = x .- x0
        
        # Use coordinates directly (already in proper units)
        xc = X[1] 
        yc = X[2] 
        zc = X[3] 
            
        # Coordinate calculations
        R = sqrt(xc^2 + yc^2 + zc^2)  # Spherical radius
        r = sqrt(xc^2 + yc^2)         # Cylindrical radius
        z = zc                        # Height
            
        # Bulge potential (using converted units)
        PHI_b = -Mb / sqrt(R^2 + b_b^2)
            
        # Disk potential (using converted units)
        PHI_d = -Md / sqrt(r^2 + (a_d + sqrt(z^2 + b_d^2))^2)
            
        # Halo potential (using converted units)
        if R < Λ
            term1 = (ν1/(γ-ν1)) * 
                        log((ν1 + (R/a_h)^(γ-ν1))/(ν1 + (Λ/a_h)^(γ-ν1)))
            term2 = (Λ/a_h)^(γ-ν1) / (ν1 + (Λ/a_h)^(γ-ν1))
            PHI_h = (Mh/a_h) * (term1 - term2)
        else
            # Outside cutoff radius
            PHI_h = -(Mh/R) * ((Λ/a_h)^γ / (1 + (Λ/a_h)^(γ-1)))
        end
            
        # Total potential (already in proper geometric units)
        PHI_total = PHI_b + PHI_d + PHI_h
            
        return PHI_total
    end
    return Φ
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
    return x->-G*M/norm(x-xo)
end #-------------------------------------------------------------------

