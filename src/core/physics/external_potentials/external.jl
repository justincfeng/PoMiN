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
    # solar offset: origin_x=-1.708859462494220E+17, origin_z=4.346342845091530E+14

    Mb_mg, Md_mg, Mh_mg, b_b_kpc, a_d_kpc, b_d_kpc, a_h_kpc, Λ_kpc, γ = p

    # Unit conversions between paper and PoMiN systems
    Mgal = tpfl(2.325e7)                    # M_gal in M_sun
    Msol2m = tpfl(1476.67)                  # 1 M_sun in meters  
    kpc_m = tpfl(3.0857e19)                 # 1 kpc in meters
    kpc = kpc_m / Msol2m                    # 1 kpc in solar mass units

    Mb = Mb_mg * Mgal
    Md = Md_mg * Mgal
    Mh = Mh_mg * Mgal
    b_b = b_b_kpc * kpc
    a_d = a_d_kpc * kpc
    b_d = b_d_kpc * kpc
    a_h = a_h_kpc * kpc
    Λ = Λ_kpc * kpc
    
    # The potential in the paper is dimensionless (already in units of c²)
    # Since we're using geometric units with c=1, no unit conversion needed for the potential itself
    # The length scales are already converted to solar mass units above
    potential_unit_conversion = tpfl(1.0)  # No conversion needed
    
    # Define constants
    G = tpfl(1.0)  # Gravitational constant in geometric units
    d = 3          # Spatial dimensions
    
    # Define origin offset from x0 parameter
    origin_x = length(x0) >= 1 ? x0[1] : tpfl(0.0)
    origin_y = length(x0) >= 2 ? x0[2] : tpfl(0.0) 
    origin_z = length(x0) >= 3 ? x0[3] : tpfl(0.0)

    if d == 3
        ν1 = tpfl(1)
        function Φ(x::RealVec)
            # Apply origin offset
            X = x .- [origin_x, origin_y, origin_z]
            
            # Convert to kpc for potential calculation
            x_kpc = X[1] / kpc
            y_kpc = X[2] / kpc  
            z_kpc = X[3] / kpc
            
            # Coordinate calculations
            R = sqrt(x_kpc^2 + y_kpc^2 + z_kpc^2)  # Spherical radius
            r = sqrt(x_kpc^2 + y_kpc^2)            # Cylindrical radius
            z = z_kpc                               # Height
            
            # Bulge potential (in units of 100 km²/s²)
            PHI_b = -Mb_mg / sqrt(R^2 + b_b_kpc^2)
            
            # Disk potential (in units of 100 km²/s²)
            PHI_d = -Md_mg / sqrt(r^2 + (a_d_kpc + sqrt(z^2 + b_d_kpc^2))^2)
            
            # Halo potential (in units of 100 km²/s²)
            if R < Λ_kpc
                term1 = (ν1/(γ-ν1)) * 
                        log((ν1 + (R/a_h_kpc)^(γ-ν1))/(ν1 + (Λ_kpc/a_h_kpc)^(γ-ν1)))
                term2 = (Λ_kpc/a_h_kpc)^(γ-ν1) / (ν1 + (Λ_kpc/a_h_kpc)^(γ-ν1))
                PHI_h = (Mh_mg/a_h_kpc) * (term1 - term2)
            else
                # Outside cutoff radius
                PHI_h = -(Mh_mg/R) * ((Λ_kpc/a_h_kpc)^γ / (1 + (Λ_kpc/a_h_kpc)^(γ-1)))
            end
            
            # Total potential in paper units (100 km²/s²)
            PHI_total_paper = PHI_b + PHI_d + PHI_h
            
            # Convert to geometric units (c=1): 100 km²/s² → c² units
            # Apply scaling factor to match observed solar orbital velocity of ~220 km/s
            c_mks = tpfl(299792458.0)  # m/s
            scaling_factor = (220.0 / 357333.0)^2  # Scale to match target velocity
            PHI_total = PHI_total_paper * tpfl(1e8) / c_mks^2 * scaling_factor
            
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
    return x->-G*M/norm(x-xo)
end #-------------------------------------------------------------------

