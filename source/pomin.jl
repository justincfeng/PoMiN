#-----------------------------------------------------------------------
#
#   PoMiN: A Hamiltonian Post-Minkowski N-Body Code
#
#       Version 2.0
#
#-----------------------------------------------------------------------

include("HamPM.jl")     # post-Minkowski Hamiltonian
include("exf.jl")       # External forces
include("gwsc.jl")      # GW strain calculator
include("integrators/epsi.jl")      # Symplectic Integrator (Tao 2016)
include("integrators/rk4i.jl")      # 4th order Runge-Kutta

module pomin

using LinearAlgebra
using CSV

import ..HPM
import ..exf
import ..epsi
import ..rk4i
import ..gwsc

const RealVec{T<:Real} = Array{T,1}

struct Particles
    m::RealVec
    q::Array{RealVec,1}
    p::Array{RealVec,1}
end

struct Parameters
    sym :: Bool
    rkl :: Tuple{Bool,Real}
    jli :: Bool
    ω :: Real
    tspan :: Tuple{Real,Real}
    iter  :: Int
    tol   :: Real
end

function Part2Z( Part::Particles )
    tpfl = typeof(Part.q[1][1])
    n = length(Part.m)
    d = length(Part.q[1])

    if n==length(Part.p) && d==length(Part.p[1])
        Z = zeros(tpfl,2*n*d)
            for i=1:n
                for j=1:d
                    Z[d*(i-1)+j]   = Part.q[i][j]
                    Z[d*(i-1+n)+j] = Part.p[i][j]
                end
            end
        return Z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end

function pominmain( Part::Particles , Param::Parameters )
    δ = abs(Param.tspan[2]-Param.tspan[1])/Param.iter

    if Param.sym
        return hsintegrator( Part2Z(Part) , Zv->dH( length(Part.q[1]) , Part.m , Zv ) , δ , Param.ω , Param.tspan , Param.iter )
    elseif Param.rkl[1]
        return hsintegrator( Part2Z(Part) , Zv->dH( length(Part.q[1]) , Part.m , Zv ) , δ , Param.rkl[2] , Param.tspan , Param.iter )
    end
    elseif Param.jli
        return hjlintegrator( Part2Z(Part) , Zv->dH( length(Part.q[1]) , Part.m , Zv ) , Param.tspan , Param.tol )
    end
end

end # module
