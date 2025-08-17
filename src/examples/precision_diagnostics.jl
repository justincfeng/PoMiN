#!/usr/bin/env julia
"""
Precision Diagnostics Script for PoMiN

This script independently tests potential precision issues:
1. ForwardDiff automatic differentiation precision
2. Mathematical function precision (sqrt, abs, etc.)
3. BigFloat vs Float64 contamination in complex calculations

Run this to isolate precision leaks from the main integration workflow.
"""

using Printf
using ForwardDiff
include("../pomin.jl")

function test_forwarddiff_precision()
    println("=== Testing ForwardDiff Precision ===")
    
    # Set BigFloat precision
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Test case: simple function that should preserve precision
    function test_func(x)
        return x[1]^2 + x[2]^2 + sqrt(x[1]^2 + x[2]^2)
    end
    
    # Create BigFloat test vector
    x_bf = [tpfl(1.23456789012345678901234567890), 
            tpfl(2.34567890123456789012345678901)]
    
    println("Input precision: $(precision(x_bf[1])) bits")
    println("Input values:")
    println("  x[1] = $(x_bf[1])")
    println("  x[2] = $(x_bf[2])")
    
    # Test function evaluation
    f_val = test_func(x_bf)
    println("\nFunction value: $(f_val)")
    println("Function value type: $(typeof(f_val))")
    println("Function value precision: $(precision(f_val)) bits")
    
    # Test ForwardDiff gradient
    grad = ForwardDiff.gradient(test_func, x_bf)
    println("\nGradient: $(grad)")
    println("Gradient type: $(typeof(grad))")
    println("Gradient element type: $(typeof(grad[1]))")
    if typeof(grad[1]) <: BigFloat
        println("Gradient precision: $(precision(grad[1])) bits")
    else
        println("WARNING: Gradient lost BigFloat precision!")
    end
    
    # Compare with manual gradient calculation
    manual_grad = [2*x_bf[1] + x_bf[1]/sqrt(x_bf[1]^2 + x_bf[2]^2),
                   2*x_bf[2] + x_bf[2]/sqrt(x_bf[1]^2 + x_bf[2]^2)]
    
    println("\nManual gradient: $(manual_grad)")
    println("Manual gradient type: $(typeof(manual_grad[1]))")
    
    # Check precision loss
    diff = abs.(grad - manual_grad)
    println("\nDifference (ForwardDiff - Manual): $(diff)")
    println("Max difference: $(maximum(diff))")
    
    if maximum(diff) > tpfl(1e-30)
        println("⚠️  WARNING: Significant difference detected - potential precision loss!")
    else
        println("✅ ForwardDiff appears to preserve precision")
    end
    
    return grad, manual_grad, diff
end

function test_mathematical_functions()
    println("\n=== Testing Mathematical Function Precision ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Test values
    test_vals = [tpfl(2.0), tpfl(4.0), tpfl(16.0), tpfl(100.0)]
    
    for val in test_vals
        println("\nTesting with value: $(val)")
        println("Input type: $(typeof(val)), precision: $(precision(val)) bits")
        
        # Test sqrt
        sqrt_result = sqrt(val)
        println("  sqrt($(val)) = $(sqrt_result)")
        println("  sqrt type: $(typeof(sqrt_result)), precision: $(precision(sqrt_result)) bits")
        
        # Test power operations
        pow2_result = val^2
        println("  $(val)^2 = $(pow2_result)")
        println("  power type: $(typeof(pow2_result)), precision: $(precision(pow2_result)) bits")
        
        # Test division
        div_result = val / tpfl(3.0)
        println("  $(val)/3 = $(div_result)")
        println("  division type: $(typeof(div_result)), precision: $(precision(div_result)) bits")
        
        # Check for precision loss by comparing with high-precision reference
        expected_sqrt = tpfl("$(sqrt(Float64(val)))")  # Convert back for comparison
        sqrt_diff = abs(sqrt_result - expected_sqrt)
        
        if sqrt_diff > tpfl(1e-25)  # Allow for some numerical difference
            println("  ⚠️  Potential sqrt precision issue: diff = $(sqrt_diff)")
        else
            println("  ✅ sqrt precision OK")
        end
    end
end

function test_hamiltonian_precision()
    println("\n=== Testing Hamiltonian Function Precision ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Create test particle system similar to MXIC_RK4
    m = [tpfl(1.0), tpfl(1.0)]  # Equal masses
    
    # Test state vector (positions and momenta)
    Z = [tpfl(-5.0), tpfl(-0.5), tpfl(0.0),   # q1
         tpfl(5.0),  tpfl(0.5),  tpfl(0.0),   # q2  
         tpfl(10.0), tpfl(0.0),  tpfl(0.0),   # p1
         tpfl(-10.0), tpfl(0.0), tpfl(0.0)]   # p2
    
    println("Test state vector type: $(typeof(Z[1]))")
    println("Test state precision: $(precision(Z[1])) bits")
    
    # Test Hamiltonian evaluation
    H_val = pomin.HamPM.H(Z, m)
    println("Hamiltonian value: $(H_val)")
    println("Hamiltonian type: $(typeof(H_val))")
    println("Hamiltonian precision: $(precision(H_val)) bits")
    
    # Test gradient calculation
    dH_val = pomin.HamPM.dH(Z, m)
    println("Gradient length: $(length(dH_val))")
    println("Gradient element type: $(typeof(dH_val[1]))")
    if typeof(dH_val[1]) <: BigFloat
        println("Gradient precision: $(precision(dH_val[1])) bits")
    else
        println("⚠️  WARNING: Gradient lost BigFloat precision!")
    end
    
    # Test individual Hamiltonian components
    println("\nTesting individual Hamiltonian terms...")
    
    # H0 (kinetic energy)
    E1 = pomin.HamPM.Enf(m[1], pomin.HamPM.psf([Z[7], Z[8], Z[9]]))
    E2 = pomin.HamPM.Enf(m[2], pomin.HamPM.psf([Z[10], Z[11], Z[12]]))
    H0 = E1 + E2
    println("H0 (kinetic): $(H0), type: $(typeof(H0))")
    
    # Check if components maintain precision
    if typeof(E1) <: BigFloat && typeof(E2) <: BigFloat
        println("✅ Kinetic energy terms maintain BigFloat precision")
    else
        println("⚠️  WARNING: Kinetic energy precision loss!")
    end
    
    return H_val, dH_val
end

function test_precision_contamination()
    println("\n=== Testing Precision Contamination Scenarios ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Test mixing BigFloat with literals
    bf_val = tpfl(2.0)
    
    println("Testing precision contamination scenarios:")
    
    # Scenario 1: Division by integer literal
    result1 = bf_val / 2
    println("BigFloat / 2 = $(result1), type: $(typeof(result1))")
    
    # Scenario 2: Division by BigFloat
    result2 = bf_val / tpfl(2)
    println("BigFloat / BigFloat(2) = $(result2), type: $(typeof(result2))")
    
    # Scenario 3: Multiplication by integer
    result3 = bf_val * 4
    println("BigFloat * 4 = $(result3), type: $(typeof(result3))")
    
    # Scenario 4: Multiplication by BigFloat
    result4 = bf_val * tpfl(4)
    println("BigFloat * BigFloat(4) = $(result4), type: $(typeof(result4))")
    
    # Scenario 5: Power with integer
    result5 = bf_val^2
    println("BigFloat^2 = $(result5), type: $(typeof(result5))")
    
    # Scenario 6: sqrt of mixed expression
    result6 = sqrt(bf_val^2 + 1)
    println("sqrt(BigFloat^2 + 1) = $(result6), type: $(typeof(result6))")
    
    result7 = sqrt(bf_val^2 + tpfl(1))
    println("sqrt(BigFloat^2 + BigFloat(1)) = $(result7), type: $(typeof(result7))")
    
    # Check for precision differences
    diff67 = abs(result6 - result7)
    println("Difference between mixed and pure BigFloat: $(diff67)")
    
    if diff67 > tpfl(1e-30)
        println("⚠️  WARNING: Precision contamination detected!")
    else
        println("✅ No significant precision contamination")
    end
end

function main()
    println("PoMiN Precision Diagnostics")
    println("=" ^ 50)
    
    # Test 1: ForwardDiff precision
    grad_fd, grad_manual, diff = test_forwarddiff_precision()
    
    # Test 2: Mathematical functions
    test_mathematical_functions()
    
    # Test 3: Hamiltonian precision
    H_val, dH_val = test_hamiltonian_precision()
    
    # Test 4: Precision contamination
    test_precision_contamination()
    
    println("\n" * ("=" ^ 50))
    println("Precision Diagnostics Complete")
    
    # Summary
    println("\n=== SUMMARY ===")
    max_diff = maximum(abs.(diff))
    if max_diff > BigFloat(1e-25)
        println("⚠️  ForwardDiff precision issue detected (max diff: $(max_diff))")
    else
        println("✅ ForwardDiff precision appears OK")
    end
    
    if typeof(dH_val[1]) <: BigFloat
        println("✅ Hamiltonian gradient maintains BigFloat precision")
    else
        println("⚠️  Hamiltonian gradient precision loss detected")
    end
    
    println("\nRun this script to isolate precision issues from the main integration.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
