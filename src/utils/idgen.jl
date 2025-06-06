#-----------------------------------------------------------------------
#
#   Initial Data Generation Functions
#
#-----------------------------------------------------------------------

"""
    setup_binary_system(m1::Real, m2::Real, r::Real, v::Real, tpfl::Type=Float64)

Set up a binary system with two particles of masses m1 and m2, separated by
distance r, with relative velocity v. Returns a Particles object.

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
    setup_circular_orbit(m1::Real, m2::Real, r::Real, tpfl::Type=Float64; G::Real=1)

Set up a binary system in a circular orbit. Returns a Particles object with
the appropriate velocity for circular motion.

Arguments:
- m1, m2: Masses of the two particles
- r: Orbital radius (separation)
- tpfl: Type for floating-point numbers (default: Float64)
- G: Gravitational constant (default: 1)

The orbital velocity is computed using v = sqrt(G(m1 + m2)/r)
"""
function setup_circular_orbit(m1::Real, m2::Real, r::Real, tpfl::Type=Float64; G::Real=1)
    # Convert inputs
    m1, m2, r, G = tpfl(m1), tpfl(m2), tpfl(r), tpfl(G)
    
    # Compute orbital velocity
    v = sqrt(G * (m1 + m2) / r)
    
    return setup_binary_system(m1, m2, r, v, tpfl)
end #-------------------------------------------------------------------

"""
    setup_massless_test_particle(m::Real, q::Vector, p::Vector, tpfl::Type=Float64)

Create a single massless test particle with specified position and momentum.

Arguments:
- m: Mass of the particle (typically 0 for photons)
- q: Position vector [x, y, z]
- p: Momentum vector [px, py, pz]
- tpfl: Type for floating-point numbers (default: Float64)
"""
function setup_massless_test_particle(m::Real, q::Vector, p::Vector, tpfl::Type=Float64)
    m = tpfl(m)
    q = convert(Vector{tpfl}, q)
    p = convert(Vector{tpfl}, p)
    
    return Particles([m], [q], [p])
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