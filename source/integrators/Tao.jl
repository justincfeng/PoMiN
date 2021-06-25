#-----------------------------------------------------------------------

module TaoMap

using LinearAlgebra

const RealVec{T<:Real} = Array{T,1}


function Zc( q::RealVec , p::RealVec)    
    tpfl=typeof(q[1])
    n = length(q)
    if n == length(p)
        Z = zeros(tpfl,2*n)        
        for i=1:n
            Z[i]     = q[i]
            Z[n+i]   = p[i]
        end
        return Z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end

function qc( Z::RealVec )    
    tpfl=typeof(Z[1])
    n4 = length(Z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        q = zeros(tpfl,n)        
        for i=1:n
            q[i]     = Z[i]
        end
        return q
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end

function pc( Z::RealVec )    
    tpfl=typeof(Z[1])
    n4 = length(Z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        p = zeros(tpfl,n)        
        for i=1:n
            p[i]     = Z[n+i]
        end
        return p
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end

function xc( Z::RealVec )    
    tpfl=typeof(Z[1])
    n4 = length(Z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        x = zeros(tpfl,n)        
        for i=1:n
            x[i]     = Z[2*n+i]
        end
        return x
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end

function yc( Z::RealVec )    
    tpfl=typeof(Z[1])
    n4 = length(Z)
    n2 = Int(n4/2)
    if iseven(n4) && iseven(n2)
        n = Int(n2/2)
        y = zeros(tpfl,n)        
        for i=1:n
            y[i]     = Z[3*n+i]
        end
        return y
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end

function ΦA( z::RealVec , dH::Function , δ::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    Z = zeros(tpfl,nz) 

    if iseven(nz) && iseven(nz/2)

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
end

function ΦB( z::RealVec , dH::Function , δ::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    Z = zeros(tpfl,nz) 

    if iseven(nz) && iseven(nz/2)

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
end

function ΦC( z::RealVec , ω::Real , δ::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    Z = zeros(tpfl,nz) 

    if iseven(nz) && iseven(nz/2)

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
end

function ΦTao2( z::RealVec , dH::Function , δ::Real , ω::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    if iseven(nz) && iseven(nz/2)

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
end

function ΦTao4( z::RealVec , dH::Function , δ::Real , ω::Real )    
    tpfl=typeof(z[1])
    nz = length(z)

    if iseven(nz) && iseven(nz/2)
        
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
end

end # TaoMap

#-----------------------------------------------------------------------
