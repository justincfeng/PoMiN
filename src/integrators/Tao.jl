#-----------------------------------------------------------------------
"""
    Zc( q::RealVec , p::RealVec)

This function returns the concatenation of the vector q and the vector 
p.
"""
function Zc( q::RealVec , p::RealVec)
    return vcat(q,p)
end #-------------------------------------------------------------------

"""
    qc( z::RealVec )

This function returns the first quarter of the `z` vector (representing 
the extended phase space coordinate vector ``z=(q,p,x,y)``), assuming 
the dimension of `z` is divisible by four.
"""
function qc( z::RealVec )
    n4 = length(z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        return z[1:n]
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------


"""
    pc( Z::RealVec )

This function returns the second quarter of the `z` vector (representing 
the extended phase space coordinate vector ``z=(q,p,x,y)``), assuming 
the dimension of `z` is divisible by four.
"""
function pc( z::RealVec )    
    n4 = length(z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        return z[n+1:2*n]
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------

"""
    xc( Z::RealVec )

This function returns the third quarter of the `z` vector (representing 
the extended phase space coordinate vector ``z=(q,p,x,y)``), assuming 
the dimension of `z` is divisible by four.
"""
function xc( z::RealVec )    
    n4 = length(z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        return z[2*n+1:3*n]
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------

"""
    yc( Z::RealVec )

This function returns the fourth quarter of the `z` vector (representing 
the extended phase space coordinate vector ``z=(q,p,x,y)``), assuming 
the dimension of `z` is divisible by four.
"""
function yc( z::RealVec )    
    n4 = length(z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        return z[3*n+1:4*n]
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------

"""
    ΦA( z::RealVec , dH::Function , δ::Real )

This function implements the function ``Φ^δ_{H_A}`` defined in the paper
arXiv:1609.02212. It takes as inputs the extended phase space vector 
`z`, a function `dH(z)` and a real-valued parameter `δ`. This function
returns a vector with the same length as the phase space vector `z`.
"""
function ΦA( z::RealVec , dH::Function , δ::Real )
    tpfl=typeof(z[1])
    nz = length(z)

    Z = zeros(tpfl,nz) 

    if iseven(nz) && iseven(Int(nz/2))

        n = Int(nz/4)

        q = qc(z)
        p = pc(z)
        x = xc(z)
        y = yc(z)

        dHH = dH(Zc(q,y))

        for i=1:n
            Z[i]     = q[i]
            Z[n+i]   = p[i] - δ*dHH[i]
            Z[2*n+i] = x[i] + δ*dHH[n+i]
            Z[3*n+i] = y[i]
        end
        return Z
    else
        print("Inputs have inconsistent dimensionality \n")
        return Z
    end
end #-------------------------------------------------------------------

"""
    ΦB( z::RealVec , dH::Function , δ::Real )

This function implements the function ``Φ^δ_{H_B}`` defined in the paper
arXiv:1609.02212. It takes as inputs the extended phase space vector 
`z`, a function `dH(z)` and a real-valued parameter `δ`. This function
returns a vector with the same length as the phase space vector `z`.
"""
function ΦB( z::RealVec , dH::Function , δ::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    Z = zeros(tpfl,nz) 

    if iseven(nz) && iseven(Int(nz/2))

        n = Int(nz/4)

        q = qc(z)
        p = pc(z)
        x = xc(z)
        y = yc(z)

        dHH = dH(Zc(x,p))

        for i=1:n
            Z[i]     = q[i] + δ*dHH[n+i]
            Z[n+i]   = p[i]
            Z[2*n+i] = x[i] 
            Z[3*n+i] = y[i] - δ*dHH[i]
        end
        return Z
    else
        print("Inputs have inconsistent dimensionality \n")
        return Z
    end
end #-------------------------------------------------------------------

"""
    ΦC( z::RealVec , ω::Real , δ::Real )

This function implements the function ``Φ^δ_{ω H_C}`` defined in the 
paper arXiv:1609.02212. It takes as inputs the extended phase space 
vector `z`, and the real-valued parameters `ω` and `δ`. This function 
returns a vector with the same length as the phase space vector `z`.
"""
function ΦC( z::RealVec , ω::Real , δ::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    Z = zeros(tpfl,nz) 

    if iseven(nz) && iseven(Int(nz/2))

        n = Int(nz/4)

        q = qc(z)
        p = pc(z)
        x = xc(z)
        y = yc(z)

        z1 = (q+x) + (q-x)*cos(2*δ*ω) + (p-y)*sin(2*δ*ω)
        z2 = (p+y) - (q-x)*sin(2*δ*ω) + (p-y)*cos(2*δ*ω)
        z3 = (q+x) - (q-x)*cos(2*δ*ω) - (p-y)*sin(2*δ*ω)
        z4 = (p+y) + (q-x)*sin(2*δ*ω) - (p-y)*cos(2*δ*ω)

        for i=1:n
            Z[i]     = z1[i]/2
            Z[n+i]   = z2[i]/2
            Z[2*n+i] = z3[i]/2
            Z[3*n+i] = z4[i]/2
        end
        
        return Z
    else
        print("Inputs have inconsistent dimensionality \n")
        return Z
    end
end #-------------------------------------------------------------------

"""
    ΦTao2( z::RealVec , dH::Function , δ::Real , ω::Real )

This function implements the map defined in Eq. (2) of arXiv:1609.02212. 
It takes as inputs the extended phase space vector `z`, a function 
`dH(z)` and the real-valued parameters `ω` and `δ`. This function 
returns a vector with the same length as the phase space vector `z`.
"""
function ΦTao2( z::RealVec , dH::Function , δ::Real , ω::Real )
    nz = length(z)
    if iseven(nz) && iseven(Int(nz/2))

        return ΦA( 
                  ΦB( 
                     ΦC(
                        ΦB( 
                           ΦA(
                              z , dH , δ/2  # Args in ΦA
                             ) , dH, δ/2    # Args in ΦB
                          ) , ω, δ    # Args in ΦC
                       ) , dH, δ/2    # Args in ΦB
                    ) , dH, δ/2    # Args in ΦA
                 )

    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------

"""
    ΦTao4( z::RealVec , dH::Function , δ::Real , ω::Real )

This function implements the map defined in Eq. (3) of arXiv:1609.02212. 
It takes as inputs the extended phase space vector `z`, a function 
`dH(z)` and the real-valued parameters `ω` and `δ`. This function 
returns a vector with the same length as the phase space vector `z`.
"""
function ΦTao4( z::RealVec , dH::Function , δ::Real , ω::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    if iseven(nz) && iseven(Int(nz/2))
        
        γ = 1/(2-2^(tpfl(1/3)))   # There is a sign error in Tao's paper

        return ΦTao2(
                     ΦTao2(
                           ΦTao2(
                                 z , dH , γ*δ , ω
                                ) , dH , (1-2*γ)*δ , ω
                          ) , dH , γ*δ , ω
                    )
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
