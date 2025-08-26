#-----------------------------------------------------------------------
#
#   GRADIENT COMPARISON: FULL HAMILTONIAN vs TEST PARTICLE HAMILTONIAN
#   
#   This script compares the gradient with respect to a single test particle
#   between the full N-body Hamiltonian H() and the test particle 
#   Hamiltonian HT(). The difference reveals the back-reaction terms
#   that are neglected in the test particle approximation.
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf, Statistics

# Use standard precision for clarity
tpfl = Float64

println("="^70)
println("GRADIENT COMPARISON: FULL vs TEST PARTICLE HAMILTONIAN")
println("="^70)

# System parameters
d = 3  # 3D space
G = 1.0  # Gravitational constant

# Background particles: Equal mass binary system
m_binary = tpfl(1.0)   # Equal masses for both binary components
m1 = m_binary          # First binary component
m2 = m_binary          # Second binary component

# Test particle: Spacecraft
m_test = tpfl(1e-10)   # Very small test mass

# Equal mass binary setup (circular orbit around center of mass)
separation = tpfl(2.0)  # Binary separation
orbital_radius = separation / 2.0  # Each star orbits at this radius

# Positions: Binary components at ±1.0 on x-axis
q1 = tpfl.([-orbital_radius, 0.0, 0.0])  # First binary component
q2 = tpfl.([+orbital_radius, 0.0, 0.0])  # Second binary component
q_test = tpfl.([0.0, 1.5, 0.2])          # Test particle offset from binary plane

# Momenta: Circular orbital motion for binary (perpendicular velocities)
# For circular orbit: v = sqrt(G*M_total/(2*r)) where M_total = 2*m_binary
v_orbital = sqrt(2.0 * m_binary / separation)  # Orbital velocity
p1 = tpfl.([0.0, -m_binary * v_orbital, 0.0])  # First component moving -y
p2 = tpfl.([0.0, +m_binary * v_orbital, 0.0])  # Second component moving +y
p_test = tpfl.([0.1, -0.05, 0.02])             # Test particle with small velocity

println("System Setup:")
println("  Background: Equal mass binary system (m1 = m2 = 1.0)")
println("  Binary separation: $(separation)")
println("  Orbital velocity: $(v_orbital)")
println("  Test particle: Spacecraft (m = 1e-10)")
println("  Positions: Binary1=[-1,0,0], Binary2=[+1,0,0], Test=[0,1.5,0.2]")
println("  Momenta: Binary1=[0,-$(v_orbital),0], Binary2=[0,+$(v_orbital),0], Test=[0.1,-0.05,0.02]")
println()

# CORRECTED COMPARISON:
# 1. Full Hamiltonian: All 3 particles interact (Binary1-Binary2-Test)
# 2. Test particle Hamiltonian: Test particle interacts with background (Binary1+Binary2)
#    but background particles still interact with each other

# Full system: [Binary1, Binary2, Test particle] - all interact
Z_full = vcat(q1, q2, q_test, p1, p2, p_test)
m_full = [m1, m2, m_test]

# Background system: [Binary1, Binary2] - they interact with each other
Z_bg = vcat(q1, q2, p1, p2)
m_bg = [m1, m2]

# Test particle system: [Test particle] - interacts with background
ZT = vcat(q_test, p_test)
mT = [m_test]

println("Corrected Setup:")
println("  Full Hamiltonian H(): 3 particles (Binary1+Binary2+Test) all interact")
println("  Test Hamiltonian HT(): Test particle + background (Binary1+Binary2 interact)")
println("  Full system Z: $(length(Z_full)) elements")
println("  Test particle ZT: $(length(ZT)) elements") 
println("  Background Z_bg: $(length(Z_bg)) elements")
println()

# Compute gradients
println("Computing gradients...")

# Full Hamiltonian gradient (all 3 particles interact)
dH_full = pomin.HamPM.dH(Z_full, m_full, d)

# Extract test particle part from full gradient
# State vector structure: [q_sun, q_jup, q_test, p_sun, p_jup, p_test]
# Test particle is 3rd particle, so indices are:
# Positions: 7-9, Momenta: 16-18
n_particles = 3
pos_indices = (2*d+1):(3*d)      # Positions of 3rd particle: indices 7-9  
mom_indices = (n_particles*d + 2*d+1):(n_particles*d + 3*d)  # Momenta: indices 16-18
test_indices = vcat(pos_indices, mom_indices)
dH_test_from_full = dH_full[test_indices]

# Test particle Hamiltonian: ONLY test particle interaction with fixed background
# The HT() function assumes background is completely fixed and unaffected
# This is the pure test particle approximation
dHT_test = pomin.HamPM.∂(ZT -> pomin.HamPM.HT(ZT, mT, Z_bg, m_bg, d), ZT)

dpF = pomin.HamPM.Z2p(n_particles,d,3,dH_full)
dqF = pomin.HamPM.Z2q(n_particles,d,3,dH_full)
dpT = pomin.HamPM.Z2p(1,d,1,dHT_test)
dqT = pomin.HamPM.Z2q(1,d,1,dHT_test)

# DEBUG: Let's check what we're actually comparing
println("DEBUG INFORMATION:")
println("="^30)

# Check the Hamiltonian values themselves
H_full_value = pomin.HamPM.H(Z_full, m_full, d)
HT_value = pomin.HamPM.HT(ZT, mT, Z_bg, m_bg, d)
println("Hamiltonian values:")
println("  H_full = $H_full_value")
println("  HT = $HT_value")
println("  Difference = $(H_full_value - HT_value)")
println()

# Check state vectors
println("State vector check:")
println("  Z_full = $Z_full")
println("  Z_bg = $Z_bg") 
println("  ZT = $ZT")
println()

# Check if background particles are identical in both cases
q1_full = pomin.HamPM.Z2q(3, d, 1, Z_full)
q2_full = pomin.HamPM.Z2q(3, d, 2, Z_full)
p1_full = pomin.HamPM.Z2p(3, d, 1, Z_full)
p2_full = pomin.HamPM.Z2p(3, d, 2, Z_full)

q1_bg = pomin.HamPM.Z2q(2, d, 1, Z_bg)
q2_bg = pomin.HamPM.Z2q(2, d, 2, Z_bg)
p1_bg = pomin.HamPM.Z2p(2, d, 1, Z_bg)
p2_bg = pomin.HamPM.Z2p(2, d, 2, Z_bg)

println("Background particle comparison:")
println("  Full H - Particle 1: q=$q1_full, p=$p1_full")
println("  Test H - Particle 1: q=$q1_bg, p=$p1_bg")
println("  Are they equal? q: $(q1_full ≈ q1_bg), p: $(p1_full ≈ p1_bg)")
println()
println("  Full H - Particle 2: q=$q2_full, p=$p2_full")
println("  Test H - Particle 2: q=$q2_bg, p=$p2_bg")
println("  Are they equal? q: $(q2_full ≈ q2_bg), p: $(p2_full ≈ p2_bg)")
println()

println("Gradient extraction:")
println("  Full H: dqF = $dqF, dpF = $dpF")
println("  Test H: dqT = $dqT, dpT = $dpT")
println()

# Compute differences using properly extracted components
dq_diff = dqF - dqT  # Position gradient differences
dp_diff = dpF - dpT  # Momentum gradient differences

println("RESULTS:")
println("="^50)

println("Full Hamiltonian gradient (test particle part):")
for i = 1:d
    @printf("  ∂H/∂q%d = %12.6e\n", i, dqF[i])
end
for i = 1:d
    @printf("  ∂H/∂p%d = %12.6e\n", i, dpF[i])
end
println()

println("Test particle Hamiltonian gradient:")
for i = 1:d
    @printf("  ∂HT/∂q%d = %12.6e\n", i, dqT[i])
end
for i = 1:d
    @printf("  ∂HT/∂p%d = %12.6e\n", i, dpT[i])
end
println()

println("Difference (Back-reaction terms):")
for i = 1:d
    @printf("  Δ(∂H/∂q%d) = %12.6e\n", i, dq_diff[i])
end
for i = 1:d
    @printf("  Δ(∂H/∂p%d) = %12.6e\n", i, dp_diff[i])
end
println()

# Compute relative differences
gradient_diff = vcat(dq_diff, dp_diff)  # Combined for compatibility with existing code
dH_test_from_full = vcat(dqF, dpF)      # Combined for compatibility
rel_diff = abs.(gradient_diff) ./ (abs.(dH_test_from_full) .+ 1e-16)

println("Relative differences:")
for i = 1:length(rel_diff)
    coord = i <= d ? "q$(i)" : "p$(i-d)"
    @printf("  |Δ(∂H/∂%s)|/|∂H/∂%s| = %12.6e\n", coord, coord, rel_diff[i])
end
println()

# Summary statistics
max_abs_diff = maximum(abs.(gradient_diff))
max_rel_diff = maximum(rel_diff)
rms_diff = sqrt(mean(gradient_diff.^2))

println("SUMMARY STATISTICS:")
println("="^30)
@printf("Maximum absolute difference: %12.6e\n", max_abs_diff)
@printf("Maximum relative difference: %12.6e\n", max_rel_diff)
@printf("RMS difference:             %12.6e\n", rms_diff)
println()

# Physical interpretation
println("ANALYSIS:")
println("="^40)
if maximum(abs.(gradient_diff)) < 1e-12
    println("• NO SIGNIFICANT DIFFERENCE found!")
    println("• This suggests the gradients are essentially identical")
    println("• Possible explanations:")
    println("  - Test particle mass is so small that back-reaction is negligible")
    println("  - The HT() function already captures the correct physics")
    println("  - There may be no true back-reaction at this level")
else
    println("• DIFFERENCE found - represents back-reaction terms:")
    println("  - How the test particle affects the background particles")
    println("  - These terms are neglected in the test particle approximation")
    println("• For m_test << m_background, differences should be small")
    println("• Larger differences indicate stronger coupling/back-reaction")
end
println()

# Let's also check what happens with a larger test mass
println("VERIFICATION WITH LARGER TEST MASS:")
println("="^50)
m_test_large = tpfl(1e-3)  # Much larger test mass
m_full_large = [m1, m2, m_test_large]
mT_large = [m_test_large]

dH_full_large = pomin.HamPM.dH(Z_full, m_full_large, d)
dH_test_from_full_large = dH_full_large[test_indices]
dHT_test_large = pomin.HamPM.∂(ZT -> pomin.HamPM.HT(ZT, mT_large, Z_bg, m_bg, d), ZT)
gradient_diff_large = dH_test_from_full_large - dHT_test_large

@printf("With m_test = %8.1e (Jupiter mass):\n", m_test_large)
@printf("  Max |Δ∇H| = %12.6e\n", maximum(abs.(gradient_diff_large)))
@printf("  Max |Δ∇H|/|∇H| = %12.6e\n", maximum(abs.(gradient_diff_large) ./ (abs.(dH_test_from_full_large) .+ 1e-16)))
println()

println("PHYSICAL INTERPRETATION:")
println("="^40)

# Test with different test particle masses
println("MASS SCALING TEST:")
println("="^30)

test_masses = [1e-15, 1e-12, 1e-10, 1e-8, 1e-6]
println("Testing how back-reaction scales with test particle mass:")
println()

for (i, m_test_scale) in enumerate(test_masses)
    m_full_scale = [m1, m2, tpfl(m_test_scale)]
    mT_scale = [tpfl(m_test_scale)]
    
    dH_full_scale = pomin.HamPM.dH(Z_full, m_full_scale, d)
    dH_test_from_full_scale = dH_full_scale[test_indices]
    
    dHT_test_scale = pomin.HamPM.∂(ZT -> pomin.HamPM.HT(ZT, mT_scale, Z_bg, m_bg, d), ZT)
    
    gradient_diff_scale = dH_test_from_full_scale - dHT_test_scale
    max_abs_diff_scale = maximum(abs.(gradient_diff_scale))
    
    @printf("m_test = %8.1e: Max |Δ∇H| = %12.6e\n", m_test_scale, max_abs_diff_scale)
end

println()
println("Expected scaling: Back-reaction ∝ m_test (linear in test mass)")
println("="^70)
