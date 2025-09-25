
#-----------------------------------------------------------------------
#   LORENTZ FACTOR
#-----------------------------------------------------------------------
"""
    γF(v,l=one(eltype(v)))

Returns the Lorentz factor γ from the velocity vector v.
"""
function γF(v,l=one(eltype(v)))
    return l/sqrt(l^2-dot(v,v))
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   γV TO v
#-----------------------------------------------------------------------
"""
    γV2v(γV)

Returns the velocity vector v from the rescaled velocity vector γV
(where γ is the Lorentz factor).
"""
function γV2v(γV,l=one(eltype(γV)))
    return γV ./ sqrt(l^2+dot(γV,γV))
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   LORENTZ TRANSFORMATION MATRIX
#-----------------------------------------------------------------------
"""
    LTM(v,l=one(eltype(v)))

Returns the Lorentz transformation matrix for a velocity vector v.
"""
function LTM(v,l=one(eltype(v)))
    β  = norm(v)
    if β < l
        γ = γF(v,l)
        Δ = γ^2 / ( l + γ )
        return [ γ       -γ*v[1]         -γ*v[2]         -γ*v[3]       ;
                -γ*v[1]  l + Δ*(v[1]^2)  Δ*v[1]*v[2]     Δ*v[1]*v[3]   ;
                -γ*v[2]  Δ*v[2]*v[1]     l + Δ*(v[2]^2)  Δ*v[2]*v[3]   ;
                -γ*v[3]  Δ*v[3]*v[1]     Δ*v[3]*v[2]     l + Δ*(v[3]^2)]
    else
        return Matrix(tpfl.(I(4)))
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   LORENTZ TRANSFORMATION FUNCTION
#-----------------------------------------------------------------------
"""
    LorentzTransform(v,l=one(eltype(v)))

Returns the Lorentz transformation function for a velocity vector v.
"""
function LorentzTransform(v,l=one(eltype(v)))
    return (X4)->LTM(v,l)*X4
end #-------------------------------------------------------------------


#-----------------------------------------------------------------------
#   THE MINKOWSKI PRODUCT (RECTANGULAR COORDINATES)
#-----------------------------------------------------------------------
"""
    η( V1::RealVec , V2::RealVec )

The function `η` takes two vectors ``V_1`` and ``V_2`` of equal
dimension, and computes their Minkowski product ``η(V_1,V_2)``, with the
assumption that the first element of the vectors corresponds to the time
component.

"""
function η( V1::RealVec , V2::RealVec )
    nv1 = length(V1)
    nv2 = length(V2)

    met = -V1[1]*V2[1]

    if nv1==nv2
        for i=2:nv1
            met += V1[i]*V2[i]
        end
        return met
    else
        print("Vectors are of a different dimension.")
        return 0*V[1]
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   THE MINKOWSKI NORM
#-----------------------------------------------------------------------
"""
    mnorm( V ) 

The function `mnorm` computes the Minkowski norm, which is
mathematically equivalent to the absolute value of the square root of
the Minkowski product ``√|η(V_1,V_2)|``. However, this function computes
the result according to the formula used in the hypot function.

"""
function mnorm( V )
    l = length(V)

    t = abs(V[1])
    x = norm(V[2:l])

    if t > x
        return abs(t - ( x/t )*( x / ( 1 + √( 1 - (x/t)^2 ) ) ))
    elseif x > t
        return abs(x - ( t/x )*( t / ( 1 + √( 1 - (t/x)^2 ) ) ))
    else
        return zero( typeof(V[1]) )
    end
end      #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   THE MINKOWSKI METRIC (RECTANGULAR COORDINATES)
#-----------------------------------------------------------------------
"""
    ημν( tpfl::DataType , dim::Int )

The function `ημν` constructs the components of the Minkowski metric.
The first argument `tpfl` specifies the floating point datatype
(typically Float64 or Double64 if one uses the DoubleFloats package) and
the second argument `dim` specifies the dimension. The default values
are `tpfl=Float64` and `dim=4`:

    julia> ημν() == ημν( Float64 , 4 )
        true

"""
function ημν( tpfl::DataType=Float64 , dim::Int=4 )   # Minkowski metric
    gη = one(tpfl)*(I(dim))

    gη[1,1] = -gη[1,1]

    return Matrix(gη)
end     #---------------------------------------------------------------