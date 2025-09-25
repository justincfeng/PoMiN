#-----------------------------------------------------------------------
#
#   Initial Data Generation Functions
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   KEPLERIAN ORBIT FUNCTIONS
#-----------------------------------------------------------------------
"""
    keplerian_orbit(a,M,e,periapsis=false)

Returns the initial conditions for a Keplerian orbit at apoapsis 
(maximum separation).

Arguments:
- a: Semi-major axis of the orbit
- M: Total mass of the system
- e: Eccentricity (0 = circular, 0 < e < 1 = elliptical)
- periapsis: Boolean indicating whether to use periapsis (default: false)
"""
function keplerian_orbit(a,M,e,periapsis=false)
    if periapsis
        r0 = a*(1-e)  # Periapsis
    else
        r0 = a*(1+e)  # Apoapsis
    end
    ve = sqrt(abs( M*(2.0/r0-1.0/a) ))
    l  = ve*r0            # Angular momentum
    ϵ = ve^2/2 - M/r0     # Total energy
    Te = 2*π*sqrt(a^3/M)
    return (r0,ve,l,ϵ,Te)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   ELEMENTARY CONSTRUCTOR FUNCTIONS
#-----------------------------------------------------------------------

"""
    setup_single_particle(m::Real, q::Vector, p::Vector, tpfl::Type=Float64)

Create a single particle with specified position and momentum.

Arguments:
- m: Mass of the particle
- q: Position vector [x, y, z]
- p: Momentum vector [px, py, pz]
- tpfl: Type for floating-point numbers (default: Float64)
"""
function setup_single_particle(m::Real, q::Vector, p::Vector, tpfl::Type=Float64)
    # Check if we're dealing with ForwardDiff dual numbers
    if eltype(p) <: Real && !(eltype(p) <: AbstractFloat && eltype(p) == tpfl)
        # For ForwardDiff compatibility: preserve the dual number type in momentum
        # but convert mass and position to target type
        m_converted = tpfl(m)
        q_converted = convert(Vector{tpfl}, q)
        # Don't convert p if it contains dual numbers - preserve for differentiation
        return Particles([m_converted], [q_converted], [p])
    else
        # Normal case: convert everything to target type
        m_converted = tpfl(m)
        q_converted = convert(Vector{tpfl}, q)
        p_converted = convert(Vector{tpfl}, p)
        return Particles([m_converted], [q_converted], [p_converted])
    end
end #-------------------------------------------------------------------

"""
    add_particle(system::Particles, m::Real, q::Vector, p::Vector)

Add a new particle to an existing particle system.

Arguments:
- system: Existing Particles object
- m: Mass of the new particle
- q: Position vector [x, y, z] of the new particle
- p: Momentum vector [px, py, pz] of the new particle

Returns a new Particles object with the additional particle.
The type of m, q, and p will be converted to match the system's type.
"""
function add_particle(system::Particles, m::Real, q::Vector, p::Vector)
    # Determine system's type from existing particles
    tpfl = eltype(system.m)
    
    # Check if we're dealing with ForwardDiff dual numbers in momentum
    if eltype(p) <: Real && !(eltype(p) <: AbstractFloat && eltype(p) == tpfl)
        # For ForwardDiff compatibility: preserve dual number type in momentum
        m_new = convert(tpfl, m)
        q_new = convert(Vector{tpfl}, q)
        # Don't convert p if it contains dual numbers
        p_new = p
    else
        # Normal case: convert everything to system type
        m_new = convert(tpfl, m)
        q_new = convert(Vector{tpfl}, q)
        p_new = convert(Vector{tpfl}, p)
    end
    
    # Create new arrays with additional particle
    masses = vcat(system.m, m_new)
    positions = vcat(system.q, [q_new])
    momenta = vcat(system.p, [p_new])
    
    return Particles(masses, positions, momenta)
end #-------------------------------------------------------------------


"""
    merge_particle_systems(systems::Vararg{Particles})

Combine multiple particle systems into a single system.
Useful for adding test particles to an existing system.
"""
function merge_particle_systems(systems::Vararg{Particles})
    masses = vcat([system.m for system in systems]...)
    positions = vcat([system.q for system in systems]...)
    momenta = vcat([system.p for system in systems]...)
    
    return Particles(masses, positions, momenta)
end #-------------------------------------------------------------------

"""
    to_com_frame(system::Particles)

Transform a particle system to the center of mass frame.
This involves:
1. Shifting positions so the center of mass is at the origin
2. Adjusting momenta so the total momentum is zero

Arguments:
- system: Particles object to transform

Returns a new Particles object in the center of mass frame.
"""
function to_com_frame(system::Particles)
    tpfl = eltype(system.m)
    N = length(system.m)  # number of particles
    d = length(system.q[1])  # number of dimensions
    
    # Calculate total mass
    M_total = sum(system.m)
    
    # Calculate center of mass position
    R_cm = zeros(tpfl, d)
    for i in 1:N
        R_cm .+= system.m[i] .* system.q[i]
    end
    R_cm ./= M_total
    
    # Calculate total momentum
    P_total = zeros(tpfl, d)
    for i in 1:N
        P_total .+= system.p[i]
    end
    
    # Create new positions relative to center of mass
    new_positions = [q .- R_cm for q in system.q]
    
    # Create new momenta in zero-momentum frame
    new_momenta = [p .- (system.m[i]/M_total) .* P_total for (i,p) in enumerate(system.p)]
    
    return Particles(copy(system.m), new_positions, new_momenta)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   SYSTEM CONSTRUCTOR FUNCTIONS
#-----------------------------------------------------------------------

"""
    setup_binary_system(m1::Real, m2::Real, r::Real, v::Real, tpfl::Type=Float64)

Set up a binary system with two particles of masses m1 and m2, separated
by distance r on the x-axis (with the center of mass at the origin), with
relative velocity v in the y-direction, in the center of mass frame.
Returns a Particles object.

Arguments:
- m1, m2: Masses of the two particles
- r: Initial separation distance
- v: Initial relative velocity magnitude
- tpfl: Type for floating-point numbers (default: Float64)

The system is set up with:
- Center of mass at origin
- Particles on x-axis
- Velocity in y-direction
- Total momentum zero (center of mass frame)
"""
function setup_binary_system(m1::Real, m2::Real, r::Real, v::Real, tpfl::Type=Float64)
    # Convert inputs to specified type
    m1, m2 = tpfl(m1), tpfl(m2)
    r, v = tpfl(r), tpfl(v)
    
    # Position vectors (along x-axis)
    x1 = -m2/(m1 + m2) * r  # COM condition
    x2 = m1/(m1 + m2) * r
    q1 = tpfl[x1, 0, 0]
    q2 = tpfl[x2, 0, 0]
    
    # Momentum vectors (along y-axis)
    p1 = tpfl[0, -m2/(m1 + m2) * v, 0] * m1  # COM frame
    p2 = tpfl[0, m1/(m1 + m2) * v, 0] * m2
    
    return Particles([m1, m2], [q1, q2], [p1, p2])
end #-------------------------------------------------------------------

"""
    setup_elliptical_orbit(m1::Real=1, m2::Real=1, a::Real=1, e::Real=0.0, 
                          tpfl::Type=Float64; G::Real=1)

Set up a binary system in an elliptical orbit at apoapsis (maximum separation).
Returns a Particles object with the appropriate position and velocity for the 
specified orbital elements in the center of mass frame.

Arguments:
- m1, m2: Masses of the two particles
- a: Semi-major axis of the orbit
- e: Eccentricity (0 = circular, 0 < e < 1 = elliptical, default: 0)
- tpfl: Type for floating-point numbers (default: Float64)
- G: Gravitational constant (default: 1)

The system is initialized at apoapsis (ν = π) where:
- Separation is maximum: r = a(1+e)
- Radial velocity is zero (turning point)
- Tangential velocity is minimum

For circular orbits (e=0), this reduces to the standard circular orbit formula.
"""
function setup_elliptical_orbit(m1::Real=1, m2::Real=1, a::Real=1, e::Real=0.0, 
                               tpfl::Type=Float64; G::Real=1)
    # Convert inputs to specified type
    m1, m2, a, e, G = tpfl(m1), tpfl(m2), tpfl(a), tpfl(e), tpfl(G)
    
    # Set initial conditions at apoapsis (ν = π, maximum separation)
    ν = tpfl(π)  # True anomaly at apoapsis
    
    # Total mass
    M = m1 + m2
    
    # Orbital radius at apoapsis
    r = a * (1 + e)  # Maximum separation
    
    # Velocity at apoapsis (purely tangential, no radial component)
    # At apoapsis: v_r = 0, v_t = sqrt(GM/a) * sqrt((1-e)/(1+e))
    v = sqrt(G * M / a) * sqrt((1 - e) / (1 + e))
    
    # Use setup_binary_system with apoapsis separation and velocity
    return setup_binary_system(m1, m2, r, v, tpfl)
end #-------------------------------------------------------------------

"""
    setup_circular_orbit(m1::Real, m2::Real, r::Real, tpfl::Type=Float64; G::Real=1)

Set up a binary system in a circular orbit. Returns a Particles object
with the appropriate velocity for circular motion.

Arguments:
- m1, m2: Masses of the two particles
- r: Orbital radius (separation)
- tpfl: Type for floating-point numbers (default: Float64)
- G: Gravitational constant (default: 1)

The orbital velocity is computed using v = sqrt(G(m1 + m2)/r)
"""
function setup_circular_orbit(m1::Real, m2::Real, r::Real, tpfl::Type=Float64; G::Real=1)
    # Convert inputs to specified type
    m1, m2 = tpfl(m1), tpfl(m2)
    r, G = tpfl(r), tpfl(G)

    μ = m1 * m2 / (m1 + m2)
    
    # Calculate total mass and circular orbital velocity for each particle
    M = m1 + m2
    v_orbital = sqrt(G * μ / r)
    
    # Position vectors (along x-axis, COM frame)
    x1 = -m2/M * r  # COM condition
    x2 = m1/M * r
    q1 = tpfl[x1, 0, 0]
    q2 = tpfl[x2, 0, 0]
    
    # Momentum vectors for circular motion (each particle has orbital velocity)
    # In COM frame, momenta are opposite but scaled by mass
    p1 = tpfl[0, -v_orbital, 0] * m1
    p2 = tpfl[0, v_orbital, 0] * m2
    
    return Particles([m1, m2], [q1, q2], [p1, p2])
end #-------------------------------------------------------------------

"""
    setup_scattering(m1::Real, m2::Real, p::Real, b::Real, d::Real; 
                     c::Real=1, tpfl::Type=Float64)

Set up initial data for a two-particle scattering problem in the center of mass frame.

Arguments:
- m1: Mass of first particle
- m2: Mass of second particle
- p: Magnitude of momentum in COM frame
- b: Impact parameter
- d: Initial separation distance
- c: Speed of light (default: 1)
- tpfl: Type for floating-point numbers (default: Float64)

The setup places:
- Particles initially separated by distance d along x-axis
- Impact parameter b split equally in y-direction
- Momenta aligned with x-axis, equal and opposite in COM frame
- For massless particles (m=0), momentum determines direction of motion

Returns a Particles object in the COM frame.
"""
function setup_scattering(m1::Real, m2::Real, p::Real, b::Real, dx::Real; c::Real=1, tpfl::Type=Float64)
    # Convert all inputs to specified type
    m1, m2 = tpfl(m1), tpfl(m2)
    p, b, dx = tpfl(p), tpfl(b), tpfl(dx)
    c = tpfl(c)  # Convert c to match system type

    # Position vectors
    q1 = tpfl[-dx/2, -b/2, 0]  # Left particle
    q2 = tpfl[dx/2, b/2, 0]    # Right particle
    
    # Momentum vectors (equal and opposite in x-direction)
    p1 = tpfl[p, 0, 0]      # Moving right
    p2 = tpfl[-p, 0, 0]     # Moving left
    
    # Create particle system
    system = Particles([m1, m2], [q1, q2], [p1, p2])
    
    # For massless particles, normalize momentum to c
    if m1 == 0 || m2 == 0
        E1 = c * sqrt((c*m1)^2 + p^2)
        E2 = c * sqrt((c*m2)^2 + p^2)
        v1 = p*c/E1
        v2 = p*c/E2
        system = Particles([m1, m2], [q1, q2], 
                          [p1 .* (v1/c), p2 .* (v2/c)])
    end
    
    return system
end #-------------------------------------------------------------------

"""
    generateUnitVectorWithinToleranceAngle(theta_tol, baseVector::RealVec, targetDistance, tpfl::Type=Float64)

Generates random unit vector that is within tolerance angle of baseVector

Arguments:
- theta_tol: Tolerance angle.  Generated unit vector will be within this angle of the baseVector
- baseVector: The vector pointing to the target
- targetDistance: Distance to the target
- tpfl: Type for floating-point numbers (default: Float64)

Returns the angle between the the unit vector and target, as well as the unit vector
"""
function generateUnitVectorWithinToleranceAngle(theta_tol, baseVector::RealVec, targetDistance, tpfl::Type=Float64)

    println("\nGenerating initial unit vector")

    diskRadius = targetDistance * tpfl(tan(deg2rad(theta_tol)))
    println("target distance = ", targetDistance * 1.47669196951425 * 6.6845871226706E-09, " AU")
    println("disk Radius = ", diskRadius * 1.47669196951425 * 6.6845871226706E-09, " AU")

    U = baseVector / norm(baseVector)

    v = nothing
    while true
        # generate vector that is linearly independent of U
        Y = [rand(tpfl), rand(tpfl), rand(tpfl)]
        # ensure Y is not collinear with U
        while abs(dot(U, Y) / (norm(U) * norm(Y))) == 1
            Y = [rand(tpfl), rand(tpfl), rand(tpfl)]
        end

        # make a vector v that is orthogonal to U by subtracting the part of Y that is parallel to U
        a = dot(Y, U) / norm(U)^2
        v = Y - a * U
        # println(stderr, "dot(v,U) = ", dot(v, U))

        # ensure v and U are orthogonal, otherwise repeat
        if dot(v, U) == 0
            break
        end
    end

    # find a vector w that's orthogonal to U and v
    w = cross(U, v)

    # scale v and w to match the disk radius
    V = v / norm(v) * diskRadius
    W = w / norm(w) * diskRadius

    # are U, V, and W all mutually orthogonal?
    # println(stderr, "U dot V = ", dot(U, V))
    # println(stderr, "U dot W = ", dot(U, W))
    # println(stderr, "V dot W = ", dot(V, W))

    # generate random point (x,y) in the unit disk
    x = nothing
    y = nothing
    while true
        x = rand(tpfl)
        y = rand(tpfl)
        if norm([x, y]) <= 1
            break
        end
    end

    # construct vector Z in the disk centered on target
    Z = x * V + y * W

    # construct init vector via vector addition of U that reaches target plus Z
    init_vec = U * targetDistance + Z

    # make init_vec a unit vector
    unit_vec = init_vec / norm(init_vec)

    # find angle between baseVector and init_vel
    cosine_theta = dot(unit_vec, baseVector) / norm(baseVector)     # init_vel already unit length
    theta_rad = acos(cosine_theta)
    theta_deg = rad2deg(theta_rad)

    println(stderr, "For unit vector ", unit_vec, " angle with base vector is ", theta_deg)

    return theta_deg, unit_vec

end

"""
    generateUnitVectorWithinToleranceAngle(theta_tol, startPt::RealVec, target::RealVec, tpfl::Type=Float64)

Generates random unit vector that points from startPt to target within a given tolerance angle

Arguments:
- theta_tol: Tolerance angle.  Generated unit vector will be within this angle of the vector pointing at target
- startPt: Starting point for unit vector that points at target
- target: Location of target
- tpfl: Type for floating-point numbers (default: Float64)

Returns the angle between the the unit vector and target, as well as the unit vector
"""
function generateUnitVectorWithinToleranceAngle(theta_tol, startPt::RealVec, target::RealVec, tpfl::Type=Float64)

    baseVector = target - startPt
    targetDistance = norm(baseVector)

    return generateUnitVectorWithinToleranceAngle(theta_tol, baseVector, targetDistance, tpfl)
end

"""
    testInitVectorsForUniformity(theta_tol, N)

Tests random initial vector generation for uniformity.  It does this by taking a target disk and dividing
    it into N rings of equal area, then counts how many initial vectors pierce each ring.  
    The total number of vectors is 1000.
    If it's uniform, the output should be a roughly equal number of vectors piercing each ring.

"""
function testInitVectorsForUniformity(theta_tol, N=6)

    startPt = Double64[-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07]  # Spacechip init pos (Moon dist (semimajor axis) from Earth in +y dir)
    target = Double64[-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13] # Proxima pos from Kervella 2017

    # test for uniformity

    theta_1 = theta_tol / sqrt(N)
    num_in_ring = zeros(N)
    for i in 1:1000
        angle = generateUnitVectorWithinToleranceAngle(theta_tol, startPt, target, Double64)[1]
        if angle <= theta_1
            num_in_ring[1] += 1
        else
            for j in 2:N
                theta_next = theta_1 * sqrt(j)
                if angle <= theta_next
                    num_in_ring[j] += 1
                    break
                end
            end
        end
    end
    println(num_in_ring)

end