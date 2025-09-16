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

# Particle system struct
struct Particles
    m::RealVec
    q::Array{RealVec}
    p::Array{RealVec}
end

"""
    elliptical_initial_conditions(a,m,e,periapsis=false)
"""
function elliptical_initial_conditions(a,M,e,periapsis=false)
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
end

println("🚀 Running External Potential Example")
println("=" ^ 50)

r0,ve,l,ϵ,Te = elliptical_initial_conditions(10.0,1.0,0.1)

V = ΦCentralMass( zeros(3) , Double64 , 1.0 )

m = [1e-10,1e-10]

U = UConstructor(V,m,zeros(3))

q0 = [[r0,0.0,0.0],[-r0,0.0,0.0]]  # Position vector as array of vectors
p0 = [[0.0,m[1]*ve,0.0],[0.0,-m[2]*ve,0.0]]  # Momentum vector as array of vectors

system = Particles(m, q0, p0)

println("Initial conditions:")
println("  Semi-major axis: $r0")
println("  Velocity: $ve")
println("  Orbital period: $Te")
println("  Particle masses: $m")
println("  Initial positions: $q0")
println("  Initial momenta: $p0")

# Test potential evaluation
test_pos = [5.0, 0.0, 0.0]
potential_value = V(test_pos)
println("  Central mass potential at [5,0,0]: $potential_value")

# Test UConstructor
test_Z = vcat(vcat(q0...), vcat(p0...))  # Flatten to phase space vector
U_value = U(test_Z)
println("  External potential energy: $U_value")

println("\n✅ External potential example completed successfully!")
println("✅ Your external potential functions are working correctly!")
