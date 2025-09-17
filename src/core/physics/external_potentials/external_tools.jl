#-----------------------------------------------------------------------
#   POTENTIAL ENERGY CONSTRUCTOR FUNCTION (NEWTONIAN)
#-----------------------------------------------------------------------

"""
    UConstructorN(PHI::Function,m::RealVec,xo::RealVec,
                              sc::Real=1.0,d::Int=3)

This function takes a potential function 'PHI(x)' and returns a new
potential energy function that is the sum of 'PHI(x)' evaluated at each
particle position 'x[i]' multiplied by the mass of the particle 'm[i]'.
It is used to construct external potentials.
"""
function UConstructorN(PHI::Function,m::RealVec,xo::RealVec,
                              sc::Real=1.0,d::Int=3)
    tpfl = typeof(m[1])
    n = length(m)

    function U(Z::RealVec)
        V = zero(eltype(Z))  # Use eltype(Z) for ForwardDiff compat.
        for i in 1:n
            x = (sc .* Z2q( n , d , i , Z )) + xo  
            V += m[i]*PHI(x)
        end
        return V
    end
    return U
end #-------------------------------------------------------------------

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
        V = zero(eltype(Z))  # Use eltype(Z) for ForwardDiff compat.
        for i in 1:n
            x = (sc .* Z2q( n , d , i , Z )) + xo  
            V += m[i]^2*PHI(x)/(Enf(m[i],psf(Z2p(n,d,i,Z))))
        end
        return V
    end
    return U
end #-------------------------------------------------------------------