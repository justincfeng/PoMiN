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

println("🔍 Comparing UConstructor vs UConstructorN")
println("=" ^ 60)

# Set up test system
r0, ve, l, ϵ, Te = let a=10.0, M=1.0, e=0.1
    r0 = a*(1+e)  # Apoapsis
    ve = sqrt(abs(M*(2.0/r0-1.0/a)))
    l = ve*r0
    ϵ = ve^2/2 - M/r0
    Te = 2*π*sqrt(a^3/M)
    (r0, ve, l, ϵ, Te)
end

# Create central mass potential
V = ΦCentralMass(zeros(3), Double64, 1.0)

# Test masses
m = [1e-10, 1e-10]

# Create both constructors
U_rel = UConstructor(V, m, zeros(3))      # Relativistic version
U_newt = UConstructorN(V, m, zeros(3))    # Newtonian version

# Test positions and momenta
q0 = [[r0, 0.0, 0.0], [-r0, 0.0, 0.0]]
p0 = [[0.0, m[1]*ve, 0.0], [0.0, -m[2]*ve, 0.0]]

# Create phase space vector
test_Z = vcat(vcat(q0...), vcat(p0...))

println("Test Configuration:")
println("  Positions: $q0")
println("  Momenta: $p0")
println("  Masses: $m")
println()

# Compare potential values
U_rel_value = U_rel(test_Z)
U_newt_value = U_newt(test_Z)

println("Potential Energy Comparison:")
println("  UConstructor (relativistic):  $(U_rel_value)")
println("  UConstructorN (Newtonian):    $(U_newt_value)")
println("  Ratio (rel/newt):             $(U_rel_value/U_newt_value)")
println()

# Let's examine the mathematical difference
println("Mathematical Analysis:")
println("UConstructorN formula: V += m[i] * PHI(x)")
println("UConstructor formula:  V += m[i]^2 * PHI(x) / (4 * Enf(m[i], psf(p[i])))")
println()

# Calculate the relativistic correction factor for each particle
for i in 1:length(m)
    pi = Z2p(length(m), 3, i, test_Z)
    psi = psf(pi)
    Eni = Enf(m[i], psi)
    correction_factor = m[i] / (4 * Eni)
    
    println("Particle $i:")
    println("  Momentum squared: $(psi)")
    println("  Energy: $(Eni)")
    println("  Relativistic correction factor: $(correction_factor)")
    println("  m[i]^2 / (4*E[i]) vs m[i]: $(correction_factor) vs $(m[i])")
end

println()
println("🎯 Key Insight:")
println("The relativistic UConstructor includes energy-dependent corrections")
println("that modify the effective coupling strength of the external potential.")
println("This explains why the orbits are no longer circular!")
