#-----------------------------------------------------------------------
#   POST-MINKOWSKIAN HAMILTONIAN
#-----------------------------------------------------------------------

# These functions compute scalar quantities that make up the Hamiltonian

"""
    psf( p::RealVec )

Momentum squared function. Returns ``p^2=\\vec{p}⋅\\vec{p}`` given
relativistic momentum ``\\vec{p}``.
"""
function psf( p::RealVec )
    return dot(p,p)
end #-------------------------------------------------------------------

"""
    Enf( m::Real, ps::Real )

Kinetic energy function. Returns ``E=\\sqrt{m^2+p^2}`` for massive
particles (m>0), and ``E=\\sqrt{p^2}=|p|`` for massless particles (m=0),
given ``m`` and ``\\vec{p}``.
"""
function Enf( m::Real, ps::Real )
    if m > zero(m)
        return sqrt(m^2+ps)
    else
        return sqrt(ps)
    end
end #-------------------------------------------------------------------

"""
    rf( qa::RealVec , qb::RealVec )

Euclidean distance function. Returns ``r_{ab}=\\sqrt{(qa-qb)⋅(qa-qb)}``,
given particle positions ``\\vec{q}_a`` and ``\\vec{q}_b``.
"""
function rf( qa::RealVec , qb::RealVec )
    return norm(qa-qb)
end #-------------------------------------------------------------------

"""
    nabf( qa::RealVec , qb::RealVec )

Separation unit vector. Returns
``\\vec{n}_{ab}=(\\vec{q}_a-\\vec{q}_b)/\\sqrt{r_{ab}}``, given particle
positions ``\\vec{q}_a`` and ``\\vec{q}_b``.
"""
function nabf( qa::RealVec , qb::RealVec )
    return (qa-qb)/rf(qa,qb)
end #-------------------------------------------------------------------

"""
    ybaf( mb::Real, qb::RealVec , qa::RealVec , pb::RealVec ,)

Returns ``y_{ba}=\\sqrt{m_b^2+(\\vec{n}_{ab}⋅\\vec{p}_b)^2}/E_b``, where
``E_b`` is Kinetic energy for particle ``b``. Inputs are particle ``b``
mass ``m_b``, particle positions ``\\vec{q}_b`` and ``\\vec{q}_a``, and
particle momenta ``\\vec{p}_b`` and ``\\vec{p}_a``.
"""
function ybaf( mb::Real, qb::RealVec , qa::RealVec , pb::RealVec )
    return sqrt(mb^2+Θabf(qb,qa,pb)^2)/Enf(mb,psf(pb))
end #-------------------------------------------------------------------

"""
    Θabf( qa::RealVec , qb::RealVec , pa::RealVec )

Returns ``Θ_{ba}=\\vec{p}_a⋅\\vec{n}_{ab}``, where ``\\vec{p}_b`` is
momentum for particle ``b``, and ``\\vec{n}_{ab}`` is unit separation
vector.
"""
function Θabf( qa::RealVec , qb::RealVec , pa::RealVec )
    return dot(pa,nabf(qb,qa))
end #-------------------------------------------------------------------

"""
    Ξabf( pa::RealVec, pb::RealVec )

Returns ``Ξ_{ba}=\\vec{p}_a⋅\\vec{p}_{b}``, where ``\\vec{p}_b`` and
``\\vec{p}_a`` are momenta for particles ``a`` and ``b``.
"""
function Ξabf( pa::RealVec, pb::RealVec )
    return dot(pa,pb)
end #-------------------------------------------------------------------

"""
    H( d::Int , m::RealVec , Z::RealVec )

Hamiltonian function. Returns the Hamiltonian for a system of particles
with masses ``m`` and positions ``Z``.
"""
function H( d::Int , m::RealVec , Z::RealVec )
    tpfl = typeof(Z[1])
    n = length(m)

    qa, qb, pa, pb  = [zeros(tpfl,d) for _ = 1:4]

    psa, psb, Ena, Enb  = [zero(tpfl) for _ = 1:4]

    rab, yba, Θab, Θba, Ξab  = [zero(tpfl) for _ = 1:5]

    H0, H1, H2, H3  = [zero(tpfl) for _ = 1:4]

    for a=1:n
        qa = Z2q(n,d,a,Z)
        pa = Z2p(n,d,a,Z)

        psa = psf(pa)
        Ena = Enf(m[a],psa)

        H0 += Ena

        for b=1:n
            if b!=a
                qb = Z2q(n,d,b,Z)
                pb = Z2p(n,d,b,Z)

                psb = psf(pb)
                Enb = Enf(m[b],psb)

                rab = rf(qa,qb)
                yba = ybaf(m[b],qb,qa,pb)
                Θab = Θabf(qa,qb,pa)
                Θba = Θabf(qb,qa,pb)
                Ξab = Ξabf(pa,pb)
                
                H1 -= Ena*Enb*( tpfl(1) + psa/(Ena^2) + psb/(Enb^2) )/(tpfl(2)*rab)
                H2 += ( tpfl(7)*Ξab - Θab*Θba )/(tpfl(4)*rab)

                if m[b] != zero(tpfl) && yba != zero(tpfl) && yba != one(tpfl)
                    H3 += ( tpfl(2)*(-tpfl(2)*(Ξab*Θba)^2 - tpfl(2)*Θab*Θba*Ξab*psb - (Θab*psb)^2 + psb*Ξab^2)/Enb^2 + 
                            tpfl(2)*(psa*Θba^2 - (Θab*Θba)^2 + tpfl(2)*Θab*Θba*Ξab - Ξab^2 + psb*Θab^2) +
                            yba*(tpfl(3)*psa*Θba^2 - (Θab*Θba)^2 + tpfl(8)*Θab*Θba*Ξab - psa*psb + tpfl(3)*psb*Θab^2)
                            ) / (tpfl(4)*Ena*Enb*rab*yba*(yba+tpfl(1))^2)
                else
                    H3 += (tpfl(3)*Enb*psb*Θab^2 - Enb*psa*(psb - tpfl(3)*Θba^2) + 
                           Enb*Θab*Θba*(-(Θab*Θba) + tpfl(8)*Ξab) + tpfl(2)*psa*psb*abs(Θba) - 
                           tpfl(2)*psb*Θab^2*abs(Θba) - tpfl(4)*Ξab^2*abs(Θba))/(tpfl(4)*Ena*rab*(Enb + abs(Θba))^2)
                end
            end
        end
    end

    return H0+H1+H2+H3

end #-------------------------------------------------------------------

"""
    dH( d::Int , m::RealVec , Z::RealVec )

Gradient of the Hamiltonian function.
"""
function dH( d::Int , m::RealVec , Z::RealVec )
    return ForwardDiff.gradient(x -> H(d, m, x), Z)
end #-------------------------------------------------------------------

"""
    Jsympl( Zarg::RealVec )

Symplectic operator. Maps output of dH to time derivative of phase space
variables.
"""
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
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------

# Right hand side of Hamilton's equations
FHE = (d,m,z)->Jsympl(dH(d,m,z))

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
