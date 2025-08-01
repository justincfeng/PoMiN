#-----------------------------------------------------------------------
#   RK4 INTEGRATOR
#-----------------------------------------------------------------------
#
#   Fourth-order Runge-Kutta integrator for Hamiltonian systems in PoMiN.
#   Uses types defined in pomin-types.jl:
#   - RealVec{T<:Real}: Type alias for Array{T} used for state vectors
#   - soln: Mutable structure for storing integration solutions
#
#   Main functions:
#   - Jsympl: Symplectic operator for Hamiltonian systems
#   - rk4map: Core RK4 integration step
#   - hrkintegrator: Main RK4 integrator with adaptive time stepping
#-----------------------------------------------------------------------

"""
    Jsympl(Zarg::RealVec)

Symplectic operator ``\\hat{J}`` for Hamiltonian systems.

# Arguments
- `Zarg::RealVec`: Phase space vector ``\\{\\vec{q}, \\vec{p}\\}``

# Returns
- `RealVec`: Result of ``\\hat{J} z`` where ``\\hat{J} = \\begin{bmatrix} 0 & I \\\\ -I & 0 \\end{bmatrix}``

# Notes
The symplectic operator transforms the gradient of the Hamiltonian into
Hamilton's equations of motion. For a phase space vector with positions
followed by momenta, this operator swaps and negates appropriately:
- Position components → momentum components
- Momentum components → negative position components

This is essential for maintaining the symplectic structure of Hamiltonian dynamics.
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

Core fourth-order Runge-Kutta integration step.

# Arguments
- `zi::RealVec`: Current phase space state vector
- `f::Function`: Function defining the system dynamics (typically `z -> Jsympl(dH(z))`)
- `δ::Real`: Time step size

# Returns
- `RealVec`: Updated phase space state after one RK4 step

# Notes
This function implements the classical fourth-order Runge-Kutta method:
```
k1 = f(zi)
k2 = f(zi + δ*k1/2)
k3 = f(zi + δ*k2/2)
k4 = f(zi + δ*k3)
zi_new = zi + δ*(k1 + 2*k2 + 2*k3 + k4)/6
```

For Hamiltonian systems, `f` is typically the composition of the symplectic
operator and the Hamiltonian gradient: `f(z) = Jsympl(∇H(z))`.
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
    hrkintegrator(d::Int, N::Int, z0::RealVec, dH::Function, δ::Real, 
                  tadapt::Function, tspan::Tuple{Real,Real}, maxit::Real)

Fourth-order Runge-Kutta integrator for Hamiltonian systems with adaptive time stepping.

# Arguments
- `d::Int`: Number of spatial dimensions (typically 3)
- `N::Int`: Number of particles in the system
- `z0::RealVec`: Initial phase space vector ``\\{\\vec{q}, \\vec{p}\\}``
- `dH::Function`: Gradient of the Hamiltonian ``\\nabla H(z)``. Takes phase space vector 
  and returns gradient vector of same dimensionality
- `δ::Real`: Initial time step size
- `tadapt::Function`: Adaptive time-stepping function. Takes parameters:
  - `δ::Real`: Current time step
  - `zi::RealVec`: Current phase space state
  - `żi::RealVec`: Time derivative of phase space state
  Returns adapted time step size
- `tspan::Tuple{Real,Real}`: Integration time span (start_time, end_time)
- `maxit::Real`: Maximum number of integration steps

# Returns
- `soln`: Solution structure containing:
  - `d::Int`: Number of dimensions
  - `N::Int`: Number of particles
  - `t::RealVec`: Time points
  - `z::Array{RealVec,1}`: Phase space trajectory
  - `zaux::Array{RealVec,1}`: Auxiliary data

# Notes
This integrator combines the classical RK4 method with adaptive time stepping
for efficient integration of post-Minkowskian Hamiltonian systems. The function
automatically constructs the symplectic dynamics `f(z) = Jsympl(dH(z))` and
integrates Hamilton's equations while adapting the time step for stability and accuracy.

Progress is printed to stderr every 1000 iterations for long integrations.
"""
function hrkintegrator(d::Int, N::Int, z0::RealVec, dH::Function, δ::Real, tadapt::Function, tspan::Tuple{Real,Real}, maxit::Real)
    tpfl = typeof(z0[1])  # tpfl = type of data stored in z0
    zi = vec(z0)

    # initialize soln data structure
    sol = soln(d,N,zeros(tpfl, 1), [zi], [zi])

    f = zx -> Jsympl(dH(zx))

    for i = 1:maxit
        if i % 1000 == 0      # print timestep number to stderr every 1000 timesteps
            println(stderr,i)   
        end
        δ = tadapt(δ,zi,f(zi))
        new_t = sol.t[i] + δ
        if new_t > tspan[2]
            break
        end
        zi = rk4map(zi, f, δ)
        sol.t = [sol.t; new_t]     # append to t
        sol.z = [sol.z; [zi] ]     # append to z
    end

    return sol
end  # End hrkintegrator
