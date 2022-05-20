#-----------------------------------------------------------------------
#   RK4 INTEGRATOR
#-----------------------------------------------------------------------

"""
    Jsympl( Zarg::RealVec )

Symplectic operator ``\\hat{J}``.  Takes vector `z` and returns ``\\hat{J} z`` where 
``\\hat{J} = \\begin{bmatrix} 0 & I \\\\ -I & 0 \\end{bmatrix}``
"""
function Jsympl(Zarg::RealVec)
    tpfl = typeof(Zarg[1])
    n2 = length(Zarg)
    Z = zeros(tpfl, n2)

    if iseven(n2)
        n = Int(round(n2 / 2, digits=0))
        for i = 1:n
            Z[i] = Zarg[n+i]
            Z[n+i] = -Zarg[i]
        end

        return Z
    else
        return Z
    end
end

"""
    rk4map(zi::RealVec, f::Function, δ::Real)

Core function in RK4 method (the RK4 equivalent of TaoMap).  Takes function `f` and initial phase space data `zi` and executes RK4 method with time step `δ`
"""
function rk4map(zi::RealVec, f::Function, δ::Real)
    tpfl = typeof(zi[1])
    n = Int(length(zi))

    k1 = zeros(tpfl, n)
    k2 = zeros(tpfl, n)
    k3 = zeros(tpfl, n)
    k4 = zeros(tpfl, n)

    k1 = f(zi)
    k2 = f(zi + δ * k1 / 2)
    k3 = f(zi + δ * k2 / 2)
    k4 = f(zi + δ * k3)

    return zi + δ * (k1 + 2 * k2 + 2 * k3 + k4) / 6
end

"""
    hrkintegrator( z0::RealVec, dH::Function , δ::Real 
                    , tadap::Function , tspan::Tuple{Real,Real} , maxit::Real )
    
4th order Runge-Kutta integrator.
Returns a struct of type "soln" (defined in pomin-types)

# Arguments
- `z0::RealVec`: Initial values phase space vector ``\\{ \\vec{q}, \\vec{p} \\}``
- `dH::Function`: Gradient of the Hamiltonian with respect to the phase space vector ``\\vec{z}``.  Takes a single parameter ``\\vec{zi}`` of the values of the phase space variables at the current time and returns a the gradient of H which is a vector with the same dimensionality as ``\\vec{zi}``
- `δ::Real`: Time step to be used by the integrator
- `tadap::Function`: Adaptive time-stepping function.  Takes a single parameter ``\\vec{zi}`` of the values of the phase space variables at the current time and returns the time step
- `tspan::Tuple{Real,Real}`: Tuple containing the start time and end time for the integration
- `maxit::Real`: Maximum number of iterations.  When this number of iterations is exceeded, integration stops.
"""
function hrkintegrator(z0::RealVec, dH::Function, δ::Real, tadap::Function, tspan::Tuple{Real,Real}, maxit::Real)
    tpfl = typeof(z0[1])

    ddof = Int(length(z0))

    tdiff = abs(tspan[2] - tspan[1])

    n = Int(round(tdiff / δ, digits=0))

    if n > abs(Int(round(maxit, digits=0)))
        n = abs(Int(round(maxit, digits=0)))
    end

    Z = soln(zeros(tpfl, n), fill(zeros(tpfl, ddof), n), fill(zeros(tpfl, ddof), n))

    zi = z0

    f = zx -> Jsympl(dH(zx))

    for i = 1:n
        zi = rk4map(zi, f, δ)
        δ = δ * tadap(zi)
        Z.t[i] = i * δ
        Z.z[i] = zi
    end

    return Z
end  # End hrkintegrator

#-----------------------------------------------------------------------
