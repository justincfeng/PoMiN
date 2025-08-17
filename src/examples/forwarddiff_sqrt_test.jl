#!/usr/bin/env julia
"""
ForwardDiff Square Root Precision Test

Test for precision issues in automatic differentiation of:
1. Sums of square root functions
2. Ratios of square root functions  
3. Complex expressions like those in the Hamiltonian

This targets potential precision loss in ForwardDiff when dealing with
sqrt functions in complex combinations, which could affect the Hamiltonian gradient.
"""

using Printf
using ForwardDiff
include("../pomin.jl")

function test_forwarddiff_sqrt_sums()
    println("=== Testing ForwardDiff on Sums of Square Roots ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Test function similar to Hamiltonian energy terms: E1 + E2
    function energy_sum(x)
        m1, m2 = tpfl(1.0), tpfl(1.0)
        p1_sq = x[1]^2 + x[2]^2 + x[3]^2
        p2_sq = x[4]^2 + x[5]^2 + x[6]^2
        E1 = sqrt(m1^2 + p1_sq)
        E2 = sqrt(m2^2 + p2_sq)
        return E1 + E2  # Sum of square roots
    end
    
    # Test vector similar to momentum values
    x_test = [tpfl(10.0), tpfl(0.0), tpfl(0.0), 
              tpfl(-10.0), tpfl(0.0), tpfl(0.0)]
    
    println("Testing energy sum: E1 + E2")
    println("Input vector: $(x_test)")
    
    # ForwardDiff gradient
    grad_fd = ForwardDiff.gradient(energy_sum, x_test)
    println("ForwardDiff gradient: $(grad_fd)")
    println("Gradient type: $(typeof(grad_fd[1]))")
    if typeof(grad_fd[1]) <: BigFloat
        println("Gradient precision: $(precision(grad_fd[1])) bits")
    end
    
    # Manual gradient calculation
    function manual_energy_gradient(x)
        m1, m2 = tpfl(1.0), tpfl(1.0)
        p1_sq = x[1]^2 + x[2]^2 + x[3]^2
        p2_sq = x[4]^2 + x[5]^2 + x[6]^2
        E1 = sqrt(m1^2 + p1_sq)
        E2 = sqrt(m2^2 + p2_sq)
        
        # dE1/dp1_i = p1_i / E1, dE2/dp2_i = p2_i / E2
        return [x[1]/E1, x[2]/E1, x[3]/E1, x[4]/E2, x[5]/E2, x[6]/E2]
    end
    
    grad_manual = manual_energy_gradient(x_test)
    println("Manual gradient: $(grad_manual)")
    
    # Compare precision
    diff = abs.(grad_fd - grad_manual)
    max_diff = maximum(diff)
    println("Max difference: $(max_diff)")
    
    if max_diff > tpfl(1e-25)
        println("⚠️  WARNING: ForwardDiff precision loss in sqrt sums!")
    else
        println("✅ ForwardDiff sqrt sum precision OK")
    end
    
    return max_diff
end

function test_forwarddiff_sqrt_ratios()
    println("\n=== Testing ForwardDiff on Ratios of Square Roots ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Test function similar to yba: sqrt(m^2 + theta^2) / sqrt(m^2 + p^2)
    function sqrt_ratio(x)
        m = tpfl(1.0)
        theta_sq = x[1]^2  # Simplified theta term
        p_sq = x[2]^2 + x[3]^2 + x[4]^2
        
        numerator = sqrt(m^2 + theta_sq)
        denominator = sqrt(m^2 + p_sq)
        return numerator / denominator
    end
    
    # Test vector
    x_test = [tpfl(0.1), tpfl(10.0), tpfl(0.0), tpfl(0.0)]
    
    println("Testing sqrt ratio: sqrt(m^2 + θ^2) / sqrt(m^2 + p^2)")
    println("Input vector: $(x_test)")
    
    # ForwardDiff gradient
    grad_fd = ForwardDiff.gradient(sqrt_ratio, x_test)
    println("ForwardDiff gradient: $(grad_fd)")
    
    # Manual gradient calculation
    function manual_ratio_gradient(x)
        m = tpfl(1.0)
        theta_sq = x[1]^2
        p_sq = x[2]^2 + x[3]^2 + x[4]^2
        
        num = sqrt(m^2 + theta_sq)
        den = sqrt(m^2 + p_sq)
        
        # d/dx[num/den] = (den * d_num/dx - num * d_den/dx) / den^2
        d_num_dx1 = x[1] / num
        d_den_dx2 = x[2] / den
        d_den_dx3 = x[3] / den
        d_den_dx4 = x[4] / den
        
        return [(den * d_num_dx1) / den^2,
                (-num * d_den_dx2) / den^2,
                (-num * d_den_dx3) / den^2,
                (-num * d_den_dx4) / den^2]
    end
    
    grad_manual = manual_ratio_gradient(x_test)
    println("Manual gradient: $(grad_manual)")
    
    # Compare precision
    diff = abs.(grad_fd - grad_manual)
    max_diff = maximum(diff)
    println("Max difference: $(max_diff)")
    
    if max_diff > tpfl(1e-25)
        println("⚠️  WARNING: ForwardDiff precision loss in sqrt ratios!")
    else
        println("✅ ForwardDiff sqrt ratio precision OK")
    end
    
    return max_diff
end

function test_forwarddiff_hamiltonian_like()
    println("\n=== Testing ForwardDiff on Hamiltonian-like Expression ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Complex function similar to H3 with multiple sqrt terms
    function hamiltonian_like(x)
        m1, m2 = tpfl(1.0), tpfl(1.0)
        
        # Extract positions and momenta
        q1 = [x[1], x[2], x[3]]
        q2 = [x[4], x[5], x[6]]
        p1 = [x[7], x[8], x[9]]
        p2 = [x[10], x[11], x[12]]
        
        # Energy terms (sqrt functions)
        E1 = sqrt(m1^2 + sum(p1.^2))
        E2 = sqrt(m2^2 + sum(p2.^2))
        
        # Distance
        r = sqrt(sum((q1 - q2).^2))
        
        # Complex expression with multiple sqrt terms
        H = E1 + E2 - (E1 * E2) / r + (E1^2 + E2^2) / (r^2)
        
        return H
    end
    
    # Test state vector similar to MXIC_RK4
    x_test = [tpfl(-5.0), tpfl(-0.5), tpfl(0.0),   # q1
              tpfl(5.0),  tpfl(0.5),  tpfl(0.0),   # q2
              tpfl(10.0), tpfl(0.0),  tpfl(0.0),   # p1
              tpfl(-10.0), tpfl(0.0), tpfl(0.0)]   # p2
    
    println("Testing Hamiltonian-like expression with multiple sqrt terms")
    
    # ForwardDiff gradient
    grad_fd = ForwardDiff.gradient(hamiltonian_like, x_test)
    println("ForwardDiff gradient length: $(length(grad_fd))")
    println("First few gradient elements: $(grad_fd[1:3])")
    println("Gradient element type: $(typeof(grad_fd[1]))")
    
    # Test precision by comparing with finite difference
    function finite_diff_gradient(f, x, h=tpfl(1e-8))
        n = length(x)
        grad = zeros(tpfl, n)
        for i in 1:n
            x_plus = copy(x)
            x_minus = copy(x)
            x_plus[i] += h
            x_minus[i] -= h
            grad[i] = (f(x_plus) - f(x_minus)) / (2*h)
        end
        return grad
    end
    
    grad_fd_approx = finite_diff_gradient(hamiltonian_like, x_test)
    println("Finite difference gradient (first few): $(grad_fd_approx[1:3])")
    
    # Compare
    diff = abs.(grad_fd - grad_fd_approx)
    max_diff = maximum(diff)
    println("Max difference vs finite difference: $(max_diff)")
    
    if max_diff > tpfl(1e-15)  # More lenient for finite difference comparison
        println("⚠️  WARNING: Significant ForwardDiff vs finite difference discrepancy!")
    else
        println("✅ ForwardDiff Hamiltonian-like precision reasonable")
    end
    
    return max_diff
end

function main()
    println("ForwardDiff Square Root Precision Test")
    println("=" ^ 50)
    
    # Test 1: Sums of square roots
    diff1 = test_forwarddiff_sqrt_sums()
    
    # Test 2: Ratios of square roots
    diff2 = test_forwarddiff_sqrt_ratios()
    
    # Test 3: Complex Hamiltonian-like expressions
    diff3 = test_forwarddiff_hamiltonian_like()
    
    println("\n" * ("=" ^ 50))
    println("ForwardDiff Square Root Precision Test Complete")
    
    println("\n=== SUMMARY ===")
    println("Max differences found:")
    println("  Sqrt sums: $(diff1)")
    println("  Sqrt ratios: $(diff2)")
    println("  Hamiltonian-like: $(diff3)")
    
    max_overall = max(diff1, diff2, diff3)
    if max_overall > BigFloat(1e-20)
        println("\n⚠️  POTENTIAL ISSUE: ForwardDiff precision loss detected!")
        println("This could contribute to the ~0.31% error plateau in MXIC_RK4.")
    else
        println("\n✅ ForwardDiff sqrt precision appears acceptable")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
