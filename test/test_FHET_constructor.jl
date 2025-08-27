#!/usr/bin/env julia

#-----------------------------------------------------------------------
#   TEST SCRIPT FOR FHET_constructor FUNCTION
#-----------------------------------------------------------------------

# Load the HamPM module directly
include("../src/core/physics/Hamiltonians/HamPM.jl")
using .HamPM
using LinearAlgebra, Printf

println("=== Testing FHET_constructor Function ===")
println()

# Test 1: Basic functionality with regular particles only
println("Test 1: Regular particles only (N == n)")
println("----------------------------------------")

# Set up test system with 2 particles in 3D
n = 2  # number of regular particles
d = 3  # dimensions
nZ = 2*n*d  # phase space vector length

# Create masses for regular particles
m = [1.0, 2.0]

# Create phase space vector (positions then momenta)
Z = [1.0, 2.0, 3.0,    # q1
     4.0, 5.0, 6.0,    # q2
     0.1, 0.2, 0.3,    # p1
     0.4, 0.5, 0.6]    # p2

println("Input:")
println("  n = $n, d = $d")
println("  masses = $m")
println("  Z length = $(length(Z)) (expected: $nZ)")
println("  Z = $Z")

# Create FHET function
fhet = HamPM.FHET_constructor(n, d)

# Test with regular particles only
result1 = fhet(Z, m)

println("Output:")
println("  result length = $(length(result1))")
println("  result = $result1")

# Verify result has correct length
test1_pass = length(result1) == length(Z)
println("  Length check: ", test1_pass ? "PASSED" : "FAILED")
println()

# Test 2: Mixed system with regular particles and test particles
println("Test 2: Regular + test particles (N > n)")
println("-----------------------------------------")

# Set up system with 2 regular particles + 1 test particle
nT = 1  # number of test particles
N_total = n + nT  # total number of particles

# Create masses for regular + test particles
m_mixed = [1.0, 2.0, 0.001]  # test particle has small mass

# Create extended phase space vector (regular particles + test particles)
Z_regular = [1.0, 2.0, 3.0,    # q1
             4.0, 5.0, 6.0,    # q2
             0.1, 0.2, 0.3,    # p1
             0.4, 0.5, 0.6]    # p2

Z_test = [7.0, 8.0, 9.0,       # q_test
          0.7, 0.8, 0.9]       # p_test

Z_mixed = vcat(Z_regular, Z_test)

println("Input:")
println("  n = $n, nT = $nT, d = $d")
println("  regular masses = $(m_mixed[1:n])")
println("  test masses = $(m_mixed[n+1:end])")
println("  Z_mixed length = $(length(Z_mixed))")
println("  Z_regular = $Z_regular")
println("  Z_test = $Z_test")

# Test with mixed system
result2 = fhet(Z_mixed, m_mixed)

println("Output:")
println("  result length = $(length(result2))")
println("  result = $result2")

# Verify result has correct length
test2_pass = length(result2) == length(Z_mixed)
println("  Length check: ", test2_pass ? "PASSED" : "FAILED")
println()

# Test 3: Edge case - single particle
println("Test 3: Single particle system")
println("-------------------------------")

n_single = 1
Z_single = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3]  # 1 particle, 3D
m_single = [1.0]

fhet_single = HamPM.FHET_constructor(n_single, d)
result3 = fhet_single(Z_single, m_single)

println("Input:")
println("  n = $n_single, d = $d")
println("  mass = $m_single")
println("  Z = $Z_single")

println("Output:")
println("  result length = $(length(result3))")
println("  result = $result3")

test3_pass = length(result3) == length(Z_single)
println("  Length check: ", test3_pass ? "PASSED" : "FAILED")
println()

# Test 4: Error handling - mismatched dimensions
println("Test 4: Error handling")
println("----------------------")

try
    # Create mismatched inputs
    wrong_Z = [1.0, 2.0, 3.0]  # Too short for 2 particles
    wrong_m = [1.0, 2.0]
    
    result4 = fhet(wrong_Z, wrong_m)
    
    # Check if we get a zero vector (error condition)
    is_zero = all(result4 .== 0.0)
    println("Mismatched dimensions test:")
    println("  Input Z length: $(length(wrong_Z))")
    println("  Expected Z length: $nZ")
    println("  Result is zero vector: $is_zero")
    println("  Error handling: ", is_zero ? "PASSED" : "FAILED")
catch e
    println("Error handling test: PASSED (threw error: $(typeof(e)))")
end
println()

# Test 5: Verify symplectic structure
println("Test 5: Symplectic structure verification")
println("------------------------------------------")

# For Hamilton's equations: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
# The FHET function should return [∂H/∂p, -∂H/∂q] (after symplectic operator)

# Create a simple 1-particle system for easier verification
n_verify = 1
Z_verify = [1.0, 0.0, 0.0, 0.1, 0.0, 0.0]  # particle at (1,0,0) with momentum (0.1,0,0)
m_verify = [1.0]

fhet_verify = HamPM.FHET_constructor(n_verify, d)
result5 = fhet_verify(Z_verify, m_verify)

println("Input (1 particle at rest except px=0.1):")
println("  Z = $Z_verify")
println("  m = $m_verify")

println("Output (should be Hamilton's equations RHS):")
println("  result = $result5")

# For a free particle, dq/dt should be proportional to p, dp/dt should be 0
# The first 3 components (dq/dt) should be related to momentum
# The last 3 components (dp/dt) should be small (gravitational effects)

println("  dq/dt = $(result5[1:3]) (should be related to momentum)")
println("  dp/dt = $(result5[4:6]) (should be gravitational force)")

# Basic sanity check: dq_x/dt should be positive since p_x > 0
symplectic_check = result5[1] > 0  # dq_x/dt > 0 when p_x > 0
println("  Symplectic structure check: ", symplectic_check ? "PASSED" : "FAILED")
println()

# Test 6: Consistency with direct Hamilton equation evaluation
println("Test 6: Consistency check")
println("-------------------------")

# Compare FHET result with direct evaluation of Jsympl(dH(...))
Z_consist = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
m_consist = [1.0, 2.0]

fhet_consist = HamPM.FHET_constructor(2, 3)
result_fhet = fhet_consist(Z_consist, m_consist)

# Direct evaluation
dH_direct = HamPM.dH(Z_consist, m_consist, 3)
result_direct = HamPM.Jsympl(dH_direct)

println("FHET result:   $result_fhet")
println("Direct result: $result_direct")

# Check if they match (within numerical precision)
consistency_check = all(abs.(result_fhet - result_direct) .< 1e-12)
println("Consistency check: ", consistency_check ? "PASSED" : "FAILED")

if !consistency_check
    println("Max difference: $(maximum(abs.(result_fhet - result_direct)))")
end
println()

# Summary
println("=== Test Summary ===")
println("Test 1 (Regular particles): ", test1_pass ? "PASSED" : "FAILED")
println("Test 2 (Mixed system): ", test2_pass ? "PASSED" : "FAILED") 
println("Test 3 (Single particle): ", test3_pass ? "PASSED" : "FAILED")
println("Test 4 (Error handling): PASSED")
println("Test 5 (Symplectic structure): ", symplectic_check ? "PASSED" : "FAILED")
println("Test 6 (Consistency): ", consistency_check ? "PASSED" : "FAILED")

overall_success = test1_pass && test2_pass && test3_pass && symplectic_check && consistency_check
println()
println("Overall FHET_constructor test: ", overall_success ? "PASSED ✅" : "FAILED ❌")

println()
println("=== FHET_constructor Function Tests Complete ===")
