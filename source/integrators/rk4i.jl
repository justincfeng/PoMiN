#-----------------------------------------------------------------------
#   RK4 INTEGRATOR
#-----------------------------------------------------------------------

function Jsympl( Zarg::RealVec )
    tpfl=typeof(Zarg[1])
    n2 = length(Zarg)
    Z = zeros(tpfl,n2)

    if iseven(n2)
        n = Int(round(n2 / 2, digits=0))
        for i=1:n
            Z[i]    = Zarg[n+i] 
            Z[n+i]  = - Zarg[i] 
        end

        return Z
    else
        return Z
    end
end

function rk4map( zi::RealVec, dH::Function , δ::Real )
    tpfl=typeof(zi[1])
    n = Int(length(zi))

    k1 = zeros(tpfl,n)
    k2 = zeros(tpfl,n)
    k3 = zeros(tpfl,n)
    k4 = zeros(tpfl,n)

    k1 = f(zi)
    k2 = f(zi+δ*k1/2)
    k3 = f(zi+δ*k2/2)
    k4 = f(zi+δ*k3)

    return zi + δ*(k1 + 2*k2 + 2*k3 + k4)/6
end

function hrkintegrator( z0::RealVec, dH::Function , δ::Real , tadap::Function , tspan::Tuple{Real,Real} , maxit::Real )
    tpfl=typeof(z0[1])

    ddof = Int(length(z0))
    
    tdiff=abs(tspan[2]-tspan[1])

    n=Int(round(tdiff/δ, digits=0))
    
    if n > abs(Int(round(maxit, digits=0)))
        n = abs(Int(round(maxit, digits=0)))
    end

    Z = soln(zeros(tpfl,n),fill(zeros(tpfl,ddof),n))

    zi = z0

    f = zx->Jsympl(dH(zi))

    for i=1:n
        zi = rk4map( zi , f , δ )
        δ  = tadap(zi)
        Z.t[i] = i*δ
        Z.z[i] = zi
    end

    return Z
end  # End hrkintegrator

#-----------------------------------------------------------------------
