using DoubleFloats
using LinearAlgebra
using Printf

# Define required types and functions
const RealVec{T<:Real} = Array{T}

# Helper functions
function Z2q(n::Int, d::Int, i::Int, Z::RealVec)
    if i > n
        error("In Z2q: particle number i=" * string(i) * 
              " exceeds number of particles n=" * string(n))
    end
    tpfl = typeof(Z[1])
    q = zeros(tpfl, d)
    for j=1:d
        index = d * (i - 1) + j
        q[j] = Z[index]
    end
    return q
end

function Z2p(n::Int, d::Int, i::Int, Z::RealVec)
    if i > n
        error("In Z2p: particle number i=" * string(i) * 
              " exceeds number of particles n=" * string(n))
    end
    tpfl = typeof(Z[1])
    p = zeros(tpfl, d)
    for j=1:d
        index = d * (i - 1) + j + d * n
        p[j] = Z[index]
    end
    return p
end

psf(p::RealVec) = dot(p, p)
Enf(m::Real, ps::Real) = m != 0 ? sqrt(m^2 + ps) : sqrt(ps)

# Include external potentials and tools
include("src/core/physics/external_potentials/external.jl")
include("src/core/physics/external_potentials/external_tools.jl")

println("🎯 Finding Circular Orbit Velocities for Relativistic Case")
println("=" ^ 60)

# Set up test system with much lower velocities
r0 = 10.0  # Fixed radius
M = 1.0    # Central mass

# Try different velocity scaling factors
velocity_factors = [0.1, 0.2, 0.3, 0.4, 0.5]

# Create central mass potential
V = ΦCentralMass(zeros(3), Double64, M)

# Test masses
m = [1e-10, 1e-10]

println("Testing different velocity scales:")
println("Radius: $r0")
println("Central mass: $M")
println("Particle masses: $m")
println()

for (i, factor) in enumerate(velocity_factors)
    # Classical circular velocity
    ve_classical = sqrt(M/r0)
    ve = factor * ve_classical
    
    # Test positions and momenta
    q0 = [[r0, 0.0, 0.0], [-r0, 0.0, 0.0]]
    p0 = [[0.0, m[1]*ve, 0.0], [0.0, -m[2]*ve, 0.0]]
    
    # Create phase space vector
    test_Z = vcat(vcat(q0...), vcat(p0...))
    
    # Calculate relativistic correction factor
    pi = Z2p(length(m), 3, 1, test_Z)
    psi = psf(pi)
    Eni = Enf(m[1], psi)
    correction_factor = m[1] / (4 * Eni)
    
    # Calculate potential energy
    U_rel = UConstructor(V, m, zeros(3))
    U_value = U_rel(test_Z)
    
    println("Factor $factor (ve = $(ve)):")
    println("  Momentum squared: $(psi)")
    println("  Energy: $(Eni)")
    println("  Relativistic correction: $(correction_factor)")
    println("  Potential energy: $(U_value)")
    println("  v/c ratio: $(ve/1.0)")  # Assuming c=1 in natural units
    println()
end

# Recommend a good velocity for circular orbits
recommended_factor = 0.2
ve_recommended = recommended_factor * sqrt(M/r0)

println("🎯 Recommendation:")
println("Use velocity factor: $recommended_factor")
println("Recommended velocity: $ve_recommended")
println("This should give more circular orbits with the relativistic UConstructor")
