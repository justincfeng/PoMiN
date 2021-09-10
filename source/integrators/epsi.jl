#-----------------------------------------------------------------------
#   EXTENDED PHASE SPACE INTEGRATOR
#-----------------------------------------------------------------------

include("Tao.jl")

function zdoubler( z::RealVec )  # This function doubles the phase space
    tpfl=typeof(z[1])
    n = length(z)
    Z = zeros(tpfl,Int(round(2*n, digits=0)))

    for i=1:n
        Z[i]     = z[i]
        Z[n+i]   = z[i]
    end    

    return Z
end  # End zdoubler

function z( Z::RealVec )  # This function extracts first half of phase space
    tpfl=typeof(Z[1])
    nf = length(Z)
    if iseven(nf)
        n = Int(round(nf/2, digits=0))
        z = zeros(tpfl,n)
        for i=1:n
            z[i] = Z[i]
        end    
        return z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end  # End zaux

function zaux( Z::RealVec )  # This function extracts second half of phase space
    tpfl=typeof(Z[1])
    nf = length(Z)
    if iseven(nf)
        n = Int(round(nf/2, digits=0))
        z = zeros(tpfl,n)
        for i=1:n
            z[i] = Z[n+i]
        end    
        return z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end  # End zaux

function hsintegrator( z0::RealVec, dH::Function , δ::Real , ω::Real , tspan::Tuple{Real,Real} , maxit::Real )
    tpfl=typeof(z0[1])

    ddof = Int(length(z0))
    
    tdiff=abs(tspan[2]-tspan[1])

    n=Int(round(tdiff/δ, digits=0))
    
    if n > abs(Int(round(maxit, digits=0)))
        n = abs(Int(round(maxit, digits=0)))
    end

    Z = soln(zeros(tpfl,n),fill(zeros(tpfl,ddof),n),fill(zeros(tpfl,ddof),n))

    zi = zdoubler(z0)

    for i=1:n
        zi = TaoMap.ΦTao4( zi , dH , δ , ω )
        Z.t[i] = i*δ
        Z.z[i] = z(zi)
        Z.zaux[i] = zaux(zi)
    end

    return Z
end  # End hintegrator

#-----------------------------------------------------------------------
