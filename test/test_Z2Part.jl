#!/usr/bin/env julia

#-----------------------------------------------------------------------
#   TEST SCRIPT FOR Z2Part FUNCTION
#-----------------------------------------------------------------------

include("../src/pomin.jl")
using .pomin
using LinearAlgebra, Printf

# Import the HamPM module to access Part2Z and Z2Part functions
import .pomin.HamPM

println("=== Testing Z2Part Function ===")
println()

# Test 1: Basic round-trip conversion (Part2Z -> Z2Part)
println("Test 1: Basic round-trip conversion")
println("-----------------------------------")

# Create test particles with Float64 precision
m_test = [1.0, 2.0, 0.5]
q_test = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
p_test = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]

original_particles = pomin.Particles(m_test, q_test, p_test)
println("Original particles created:")
println("  Masses: ", original_particles.m)
println("  Positions: ", original_particles.q)
println("  Momenta: ", original_particles.p)
println()

# Convert to phase space vector
Z = HamPM.Part2Z(original_particles)
println("Phase space vector Z:")
println("  Length: ", length(Z))
println("  Values: ", Z)
println()

# Convert back to particles
n = length(m_test)
d = 3
reconstructed_particles = HamPM.Z2Part(Z, m_test, n, d)
println("Reconstructed particles:")
println("  Masses: ", reconstructed_particles.m)
println("  Positions: ", reconstructed_particles.q)
println("  Momenta: ", reconstructed_particles.p)
println()

# Check if reconstruction is exact
mass_match = original_particles.m == reconstructed_particles.m
pos_match = all([original_particles.q[i] ≈ reconstructed_particles.q[i] for i in 1:n])
mom_match = all([original_particles.p[i] ≈ reconstructed_particles.p[i] for i in 1:n])

println("Verification:")
println("  Masses match: ", mass_match)
println("  Positions match: ", pos_match)
println("  Momenta match: ", mom_match)
println("  Overall success: ", mass_match && pos_match && mom_match)
println()

# Test 2: Different precision types
println("Test 2: Different precision types")
println("---------------------------------")

using DoubleFloats

# Test with DoubleFloats
m_double = Double64.([1.0, 0.5])
q_double = [Double64.([1.1, 2.2, 3.3]), Double64.([4.4, 5.5, 6.6])]
p_double = [Double64.([0.11, 0.22, 0.33]), Double64.([0.44, 0.55, 0.66])]

particles_double = pomin.Particles(m_double, q_double, p_double)
Z_double = HamPM.Part2Z(particles_double)
reconstructed_double = HamPM.Z2Part(Z_double, m_double, 2, 3)

double_success = (particles_double.m == reconstructed_double.m) &&
                 all([particles_double.q[i] ≈ reconstructed_double.q[i] for i in 1:2]) &&
                 all([particles_double.p[i] ≈ reconstructed_double.p[i] for i in 1:2])

println("DoubleFloats precision test: ", double_success ? "PASSED" : "FAILED")
println()

# Test 3: Edge cases
println("Test 3: Edge cases")
println("------------------")

# Single particle
m_single = [1.0]
q_single = [[1.0, 2.0, 3.0]]
p_single = [[0.1, 0.2, 0.3]]

particles_single = pomin.Particles(m_single, q_single, p_single)
Z_single = HamPM.Part2Z(particles_single)
reconstructed_single = HamPM.Z2Part(Z_single, m_single, 1, 3)

single_success = (particles_single.m == reconstructed_single.m) &&
                 (particles_single.q[1] ≈ reconstructed_single.q[1]) &&
                 (particles_single.p[1] ≈ reconstructed_single.p[1])

println("Single particle test: ", single_success ? "PASSED" : "FAILED")

# Two dimensions
m_2d = [1.0, 2.0]
q_2d = [[1.0, 2.0], [3.0, 4.0]]
p_2d = [[0.1, 0.2], [0.3, 0.4]]

particles_2d = pomin.Particles(m_2d, q_2d, p_2d)
Z_2d = HamPM.Part2Z(particles_2d)
reconstructed_2d = HamPM.Z2Part(Z_2d, m_2d, 2, 2)

twod_success = (particles_2d.m == reconstructed_2d.m) &&
               all([particles_2d.q[i] ≈ reconstructed_2d.q[i] for i in 1:2]) &&
               all([particles_2d.p[i] ≈ reconstructed_2d.p[i] for i in 1:2])

println("2D particles test: ", twod_success ? "PASSED" : "FAILED")
println()

# Test 4: Error handling
println("Test 4: Error handling")
println("----------------------")

try
    # Test with mismatched dimensions
    wrong_Z = [1.0, 2.0, 3.0]  # Too short
    wrong_m = [1.0, 2.0]
    result = HamPM.Z2Part(wrong_Z, wrong_m, 2, 3)
    println("Error handling test: FAILED (should have thrown error)")
catch e
    println("Error handling test: PASSED (correctly threw error: ", typeof(e), ")")
end
println()

# Test 5: Complex case with n < nF (if implemented)
println("Test 5: Complex case (n < nF)")
println("------------------------------")

try
    # Create a system with more masses than we want to extract
    m_complex = [1.0, 2.0, 3.0, 4.0]  # 4 masses
    # But we'll only extract first 2 particles
    n_extract = 2
    
    # Create phase space vector for all 4 particles
    q_complex = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
    p_complex = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9], [1.0, 1.1, 1.2]]
    
    particles_complex = pomin.Particles(m_complex, q_complex, p_complex)
    Z_complex = HamPM.Part2Z(particles_complex)
    
    # Try to extract only first 2 particles
    result = HamPM.Z2Part(Z_complex, m_complex, n_extract, 3)
    
    if isa(result, Tuple)
        part1, part2 = result
        println("Complex case returned tuple:")
        println("  First system - Masses: ", part1.m, ", Particles: ", length(part1.m))
        println("  Second system - Masses: ", part2.m, ", Particles: ", length(part2.m))
        
        # Verify first system matches first 2 particles
        first_match = (part1.m == m_complex[1:2]) &&
                      all([part1.q[i] ≈ q_complex[i] for i in 1:2]) &&
                      all([part1.p[i] ≈ p_complex[i] for i in 1:2])
        
        # Verify second system matches last 2 particles  
        second_match = (part2.m == m_complex[3:4]) &&
                       all([part2.q[i] ≈ q_complex[i+2] for i in 1:2]) &&
                       all([part2.p[i] ≈ p_complex[i+2] for i in 1:2])
        
        println("  First system match: ", first_match ? "PASSED" : "FAILED")
        println("  Second system match: ", second_match ? "PASSED" : "FAILED")
    else
        println("Complex case returned single system: ", typeof(result))
    end
catch e
    println("Complex case test: ERROR (", typeof(e), ": ", e, ")")
end
println()

println("=== Z2Part Function Tests Complete ===")
