#-----------------------------------------------------------------------
#   ADAPTIVE TIMESTEPPING
#-----------------------------------------------------------------------

module tadap

using LinearAlgebra

include("../core/pomin-types.jl")
include("../physics/Hamiltonians/idxer.jl")

#-----------------------------------------------------------------------

"""
    rf(qa::RealVec, qb::RealVec)

Compute the relativistic distance between two position vectors in the post-Minkowskian framework.
"""
function rf(qa::RealVec, qb::RealVec)
    return norm(qa - qb)
end

#-----------------------------------------------------------------------
#   SIMPLE ADAPTIVE TIMESTEPPING
#-----------------------------------------------------------------------
"""
    tcour( dt::Real , Z::RealVec , Zdot::RealVec , C = 0.001 , d=3 )

This function implements a simple adaptive timestepping function 
inspired by the Courant-Friedrichs Condition. 

Returns a new timestep that satisfies the CFL condition.
If the current timestep already satisfies CFL condition, the timestep is returned unchanged.

# Arguments
- `dt::Real` is the current timestep
- `Z::RealVec` is the current state of the system
- `Zdot::RealVec` is the time derivative of Z
- `C` is the Courant factor
- `d` is the number of dimensions

"""
function tcour( dt::Real , Z::RealVec , Zdot::RealVec , C = 0.001 , d=3 )
    tpfl = typeof(Z[1])
    n2 = length(Z)
    rs = tpfl(1)
    vs = tpfl(1)

    ndof = Int(round(n2 / 2, digits=0))
    n    = Int(round(ndof/d, digits=0))

    if iseven(n2) && n*d==ndof
        # calculate separation and relative velocity between first two particles
        # and store these as starting minimum values
        qa  = Z2q(n,d,1,Z)
        qb  = Z2q(n,d,2,Z)
        va  = Z2q(n,d,1,Zdot)
        vb  = Z2q(n,d,2,Zdot)
        vab = norm(va-vb)
        rab = rf(qa,qb)
        rs  = rab  # smallest separation among particle pairs
        vs  = vab       # relative velocity for particle pair with smallest separation
        for a=1:n
            qa = Z2q(n,d,a,Z)    
            va = Z2q(n,d,a,Zdot)
            for b=2:n
                if b!=a
                    qb = Z2q(n,d,b,Z)
                    vb = Z2q(n,d,b,Zdot)
                    rab = rf(qa,qb)
                    vab = norm(va-vb)
                    if rab<=rs  
                        vs = vab
                        rs = rab
                    end
                end
            end
        end
        Δt = tpfl(C) * rs/vs
        # return \Delta t if it is smaller than the maximum timestep, otherwise return maximum timestep

        # println("time step = " * string(Δt))
        return Δt
    else
        error("in tcour delta t is not changing bc ( iseven(n2) && n*d==ndof ) returned false")
        return dt
    end
end  # End tcour

function dt_lin_int(dt::Real , Z::RealVec , Zdot::RealVec , start_dist::Real , min_dist::Real , max_dt::Real , min_dt::Real, d = 3)
"""
    This function implements a simple adaptive timestepping function by linear interpolation.
    When distance between first two bodies is equal to or greater than start_dist, it uses the max_dt
    When distance is equal to or less than min_dist, it uses the min_dt
    Between the two, it interpolates linearly
"""
    # find distance between first two bodies

    tpfl = typeof(Z[1])
    n2 = length(Z)

    ndof = Int(round(n2 / 2, digits=0))
    n    = Int(round(ndof/d, digits=0))
    if iseven(n2) && n*d==ndof
        qa  = Z2q(n,d,1,Z)
        qb  = Z2q(n,d,2,Z)
        rab = rf(qa,qb)
    else
        error("in dt_lin_int delta t is not changing bc ( iseven(n2) && n*d==ndof ) returned false")
        return dt
    end

    if rab >= start_dist
        dt = max_dt
    elseif rab <= min_dist
        dt = min_dt
    else
        # interpolate linearly
        dt = abs(start_dist - rab) / abs(start_dist - min_dist) * min_dt + abs(rab - min_dist) / abs(start_dist - min_dist) * max_dt
    end

    # double-check that we haven't gone outside range of [min_dt, max_dt]
    if dt > max_dt
        dt = max_dt
    end
    if dt <= min_dt
        dt = min_dt
    end

    return dt

end

end # end of HamPM
