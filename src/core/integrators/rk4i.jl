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
#   - rk4map: Core RK4 integration step
#   - hrkintegrator: Main RK4 integrator with adaptive time stepping
#
#   Uses Jsympl from HamTools.jl for symplectic operations.
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   RK4 MAP
#-----------------------------------------------------------------------
"""
    rk4map(zi::RealVec, f::Function, δ::Real)

Core fourth-order Runge-Kutta integration step.

# Arguments
- `zi::RealVec`: Current phase space state vector
- `f::Function`: Function defining the system dynamics
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

For Hamiltonian systems, `f` is typically the composition of the
symplectic operator and the Hamiltonian gradient: `f(z) =
Jsympl(∇H(z))`.
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
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   SIMPLE ADAPTIVE TIMESTEPPING
#-----------------------------------------------------------------------
"""
    tcour( dt::Real , Z::RealVec , Zdot::RealVec , C = 0.001 , d=3 )

This function implements a simple adaptive timestepping function 
inspired by the Courant–Friedrichs–Lewy (CFL) Condition. 

Returns a new timestep that satisfies the CFL condition. If the current
timestep already satisfies CFL condition, the timestep is returned
unchanged.

Used with the `hrkintegrator` function.

# Arguments
- `dt::Real` is the current timestep
- `Z::RealVec` is the current state of the system
- `Zdot::RealVec` is the time derivative of Z
- `C` is the Courant factor
- `d` is the number of dimensions

"""
function tcour( dt::Real , Z::RealVec , Zdot::RealVec , C = 0.001 , d=3)
    tpfl = typeof(Z[1])
    n2 = length(Z)
    rs = tpfl(1)
    vs = tpfl(1)

    ν0, ν1, ν2, ν3, ν4, ν5, ν6, ν7, ν8, ν9 = tpnum(tpfl)

    ndof = Int(round(n2 / 2, digits=0))
    n    = Int(round(ndof/d, digits=0))

    if iseven(n2) && n*d==ndof
        # calculate separation and relative velocity between first two
        # particles and store these as starting minimum values
        qa  = Z2q(n,d,1,Z)
        qb  = Z2q(n,d,2,Z)
        va  = Z2q(n,d,1,Zdot)
        vb  = Z2q(n,d,2,Zdot)
        vab = norm(va-vb)
        rab = norm(qa - qb)
        rs  = rab  # smallest separation among particle pairs
        vs  = vab  # relative velocity for closest particle pair
        for a=1:n
            qa = Z2q(n,d,a,Z)    
            va = Z2q(n,d,a,Zdot)
            for b=2:n
                if b!=a
                    qb = Z2q(n,d,b,Z)
                    vb = Z2q(n,d,b,Zdot)
                    rab = norm(qa - qb)
                    vab = norm(va-vb)
                    if rab<=rs  
                        vs = vab
                        rs = rab
                    end
                end
            end
        end
        Δt = tpfl(C) * rs/vs
        # return \Delta t if it is smaller than the maximum timestep,
        # otherwise return maximum timestep
        return Δt
    else
        error("In tcour dt unchanged bc ( iseven(n2) && n*d==ndof ) returned false")
        return dt
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   FIXED TIMESTEPPING
#-----------------------------------------------------------------------
"""
    tnone( dt::Real , Z::RealVec , Zdot::RealVec , C = 0.001 , d=3 )

This function implements a fixed timestep.
"""
function tnone(dt,z,zdot) 
    return dt
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   RK4 INTEGRATOR
#-----------------------------------------------------------------------
"""
    hrkintegrator(d::Int, N::Int, z0::RealVec, dH::Function, δ::Real, 
                       tspan::Tuple{Real,Real}, maxit::Real, 
                       tadapt::Function=tnone, Nrec::Int=100)

Fourth-order Runge-Kutta integrator for Hamiltonian systems with
adaptive time stepping.

# Arguments
- `d::Int`: Number of spatial dimensions (typically 3)
- `N::Int`: Number of particles in the system
- `z0::RealVec`: Initial phase space vector ``\\{\\vec{q}, \\vec{p}\\}``
- `dH::Function`: Gradient of the Hamiltonian ``\\nabla H(z)``. Takes
  phase space vector and returns gradient vector of same dimensionality
- `δ::Real`: Initial time step size
- `tadapt::Function`: Adaptive time-stepping function. Takes parameters:
  - `δ::Real`: Current time step
  - `zi::RealVec`: Current phase space state
  - `żi::RealVec`: Time derivative of phase space state
  Returns adapted time step size
- `tspan::Tuple{Real,Real}`: Time span (start_time, end_time)
- `maxit::Real`: Maximum number of integration steps
- `Nrec::Int`: Recording frequency (Nrec=-1 saves only last point,
  Nrec=-N saves last N points)

# Returns
- `soln`: Solution structure containing:
  - `d::Int`: Number of dimensions
  - `N::Int`: Number of particles
  - `t::RealVec`: Time points
  - `z::Array{RealVec,1}`: Phase space trajectory
  - `zaux::Array{RealVec,1}`: Auxiliary data

# Notes
This integrator combines the classical RK4 method with adaptive time
stepping for efficient integration of post-Minkowskian Hamiltonian
systems. The function automatically constructs the symplectic dynamics
`f(z) = Jsympl(dH(z))` and integrates Hamilton's equations while
adapting the time step for stability and accuracy.

Progress is printed to stderr every 1000 iterations for long
integrations.
"""
function hrkintegrator(d::Int, N::Int, z0::RealVec, dH::Function, δ::Real, 
                       tspan::Tuple{Real,Real}, maxit::Real, 
                       tadapt::Function=tnone, Nrec::Int=100)
    tpfl = typeof(z0[1])  # tpfl = type of data stored in z0
    zi = vec(z0)

    # initialize soln data structure
    sol = soln(d,N,[tpfl(tspan[1])], [zi], [zi])

    f = zx -> Jsympl(dH(zx))
    
    # Handle different Nrec cases
    if Nrec == -1
        # Save only the last point
        temp_t = [sol.t[1]]
        temp_z = [zi]
        
        for i = 1:maxit
            if i % 1000 == 0      # print timestep to stderr
                println(stderr,i)   
            end
            δ = tadapt(δ,zi,f(zi))
            new_t = temp_t[end] + δ
            if new_t > tspan[2]
                break
            end
            zi = rk4map(zi, f, δ)
            # Only keep the latest values
            temp_t = [new_t]
            temp_z = [zi]
        end
        
        # Set final solution with only initial and final states
        sol.t = [sol.t[1], temp_t[1]]
        sol.z = [sol.z[1], temp_z[1]]
        
    elseif Nrec < 0
        # Save last |Nrec| points - collect all then truncate
        all_t = [sol.t[1]]
        all_z = [zi]
        
        for i = 1:maxit
            if i % 1000 == 0      # print timestep number to stderr
                println(stderr,i)   
            end
            δ = tadapt(δ,zi,f(zi))
            new_t = all_t[end] + δ
            if new_t > tspan[2]
                break
            end
            zi = rk4map(zi, f, δ)
            all_t = [all_t; new_t]
            all_z = [all_z; [zi]]
        end
        
        # Keep only the last |Nrec| points
        n_keep = min(abs(Nrec), length(all_t))
        sol.t = all_t[end-n_keep+1:end]
        sol.z = all_z[end-n_keep+1:end]
        
    else
        # Default behavior: save every Nrec steps (or all if Nrec=1)
        step_count = 0
        current_t = sol.t[end]
        
        for i = 1:maxit
            δ = tadapt(δ,zi,f(zi))
            new_t = current_t + δ
            if new_t > tspan[2]
                break
            end
            zi = rk4map(zi, f, δ)
            current_t = new_t
            step_count += 1
            
            # Save every Nrec steps (or every step if Nrec=1)
            if step_count % max(Nrec, 1) == 0
                println(stderr,i) 
                sol.t = [sol.t; current_t]     # append to t
                sol.z = [sol.z; [zi] ]     # append to z
            end
        end
        
        # Always save the final state if we didn't just save it
        if step_count % max(Nrec, 1) != 0
            sol.t = [sol.t; current_t]
            sol.z = [sol.z; [zi]]
        end
    end

    return sol
end  # End hrkintegrator
