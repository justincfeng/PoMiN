#-----------------------------------------------------------------------
#
#   Initial Data Generation Functions
#
#-----------------------------------------------------------------------

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
    m = tpfl(m)
    q = convert(Vector{tpfl}, q)
    p = convert(Vector{tpfl}, p)
    
    return Particles([m], [q], [p])
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
    
    # Convert inputs to match system type
    m_new = convert(tpfl, m)
    q_new = convert(Vector{tpfl}, q)
    p_new = convert(Vector{tpfl}, p)
    
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
