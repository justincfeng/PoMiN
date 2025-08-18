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
                origin_x=tpfl(-1.708859462494220E+17), origin_y=tpfl(0), 
                origin_z =tpfl(4.346342845091530E+14), r_sun=tpfl(8.4), 
                Mb=tpfl(409), Md=tpfl(2856), Mh=tpfl(1018), 
                b_b=tpfl(0.23), a_d=tpfl(4.22), b_d=tpfl(0.292), 
                a_h=tpfl(2.562), Λ=tpfl(100), γ=tpfl(2.02),d::Int=3 )
    # (origin_x, origin_y, origin_z) is the location (in units of Msol)
    # in the galactocentric frame where the simulation's origin is
    # placed it defaults to the location of the Sun as given by astropy
    # and 2019 data from arXiv:1904.05721
    
    # Milky Way Model I (Irrgang et al., arXiv:1211.4353):
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

    Mgal = tpfl(2.325e7)            # mass of the Milky Way in Msol
    kpc  = tpfl(5.224206385e15)     # kpc in units of Msol

    if d == 3
        ν1 = tpfl(1)
        function Φ(x::RealVec)
                X  = x .+ [origin_x, origin_y, origin_z]  # Remove type conversion for ForwardDiff compatibility
                R     = norm(X)*kpc
                r     = sqrt(X[1]^2 + X[2]^2)*kpc
                z     = X[3]*kpc
                PHI_b = -Mb/sqrt(R^2+b_b^2)
                PHI_d = -Md/sqrt(r^2 + (a_d + sqrt(z^2 + b_d^2))^2)
                PHI_h = (Mh/a_h) * 
                    ( 
                     (ν1/(γ-ν1))*
                     log((ν1 + (R/a_h)^(γ-ν1))/
                     (ν1 + (Λ/a_h)^(γ-ν1))) - 
                     (Λ/a_h)^(γ-ν1) / (ν1 + (Λ/a_h)^(γ-ν1)) 
                    )
    
            return PHI_b + PHI_d + PHI_h
        end
        return Φ
    else 
        return x->zero(eltype(x))  # Use eltype(x) instead of tpfl for type compatibility
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

