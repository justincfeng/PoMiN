#-----------------------------------------------------------------------
#   JULIA INTEGRATORS
#-----------------------------------------------------------------------

module jli      # jli module

using LinearAlgebra
using OrdinaryDiffEq

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

function hjlintegrator( z0::RealVec , dH::Function , tspan::Tuple{Real,Real} , tol::Real )
    function f(zx,p,t)
        return Jsympl(dH(zx))
    end

    prob = ODEProblem(f,z0,tspan,abstol=tol,reltol=tol)

    S = solve(prob,Tsit5())

    return (soln(S.t,S.u),S)
end

end  # End scope of module jli
#-----------------------------------------------------------------------
