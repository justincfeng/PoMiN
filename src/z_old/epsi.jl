#-----------------------------------------------------------------------
#   EXTENDED PHASE SPACE INTEGRATOR
#-----------------------------------------------------------------------

include("Tao.jl")

"""
    zdoubler( z::RealVec )

Doubles the phase space.  Takes a phase space vector ``\\vec{z}`` and returns a new vector ``\\{\\vec{z},\\vec{z}\\}``
"""
function zdoubler(z::RealVec)  # This function doubles the phase space
    tpfl = typeof(z[1])
    n = length(z)
    Z = zeros(tpfl, Int(round(2 * n, digits=0)))

    for i = 1:n
        Z[i] = z[i]
        Z[n+i] = z[i]
    end

    return Z
end  # End zdoubler

"""
    z( Z::RealVec )

Extracts first half of even-dimensional vector ``\\vec{Z}``
"""
function z(Z::RealVec)  # This function extracts first half of phase space
    tpfl = typeof(Z[1])
    nf = length(Z)
    if iseven(nf)
        n = Int(round(nf / 2, digits=0))
        z = zeros(tpfl, n)
        for i = 1:n
            z[i] = Z[i]
        end
        return z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end  # End zaux

"""
    zaux( Z::RealVec )

Extracts second half of an even-dimensional vector ``\\vec{Z}``
"""
function zaux(Z::RealVec)  # This function extracts second half of phase space
    tpfl = typeof(Z[1])
    nf = length(Z)
    if iseven(nf)
        n = Int(round(nf / 2, digits=0))
        z = zeros(tpfl, n)
        for i = 1:n
            z[i] = Z[n+i]
        end
        return z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end  # End zaux

"""
    hsintegrator( z0::RealVec, dH::Function , δ::Real , ω::Real , tspan::Tuple{Real,Real} , maxit::Real )

Symplectic integrator of Tao.
Returns a struct of type "soln" (defined in pomin-types)

# Arguments
- `d::Int` is the number of dimensions
- `N::Int` is the number of particles
- `z0::RealVec`: Initial values of phase space vector ``\\{ \\vec{q}, \\vec{p} \\}``
- `dH::Function`: Gradient of the Hamiltonian with respect to the phase space vector ``\\vec{z}``.  Takes a single parameter ``\\vec{zi}`` of the values of the phase space variables at the current time and returns a the gradient of H which is a vector with the same dimensionality as ``\\vec{zi}``
- `δ::Real`: Time step to be used by the integrator
- `ω::Real`: Parameter needed for Tao's method of integration.  See arXiv:1609.02212
- `tspan::Tuple{Real,Real}`: Tuple containing the start time and end time for the integration
- `maxit::Real`: Maximum number of iterations.  When this number of iterations is exceeded, integration stops.
"""
function hsintegrator(d::Int, N::Int, z0::RealVec, dH::Function, δ::Real, ω::Real, tspan::Tuple{Real,Real}, maxit::Real)
    tpfl = typeof(z0[1])

    ddof = Int(length(z0))

    tdiff = abs(tspan[2] - tspan[1])

    n = Int(round(tdiff / δ, digits=0))

    if n > abs(Int(round(maxit, digits=0)))
        n = abs(Int(round(maxit, digits=0)))
    end

    #Z = soln(zeros(tpfl, n), fill(zeros(tpfl, ddof), n), fill(zeros(tpfl, ddof), n))
    Z = soln(d,N,zeros(tpfl, n), fill(zeros(tpfl, ddof), n), fill(zeros(tpfl, ddof), n))

    zi = zdoubler(vec(z0))

    for i = 1:n
        zi = ΦTao4(zi, dH, δ, ω)
        Z.t[i] = i * δ
        Z.z[i] = z(zi)
        Z.zaux[i] = zaux(zi)
    end

    return Z
end  # End hintegrator

#-----------------------------------------------------------------------
