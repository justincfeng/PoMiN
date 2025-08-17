#!/usr/bin/env julia
"""
Square Root Cancellation Diagnostics

Test for numerical cancellation issues in expressions like:
1. yba + 1 ≈ 0 (catastrophic cancellation)
2. (yba + 1)^2 in denominators (division by very small numbers)
3. sqrt operations with potential precision loss

This targets the specific issue in H3 Hamiltonian: (yba+ν1)^2
"""

using Printf
include("../pomin.jl")

function test_yba_cancellation()
    println("=== Testing yba Cancellation Issues ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    ν0, ν1, ν2, ν3, ν4, ν5, ν6, ν7, ν8, ν9 = pomin.tpnum(tpfl)
    
    # Test parameters similar to MXIC_RK4 setup
    m = [tpfl(1.0), tpfl(1.0)]  # Equal masses
    
    # Test different impact parameters to see when yba ≈ -1
    impact_params = [tpfl(10.0), tpfl(100.0), tpfl(1000.0), tpfl(10000.0), tpfl(100000.0)]
    
    for b in impact_params
        println("\n--- Testing impact parameter b = $(b) ---")
        
        # Setup scattering configuration
        p = tpfl(10.0)  # momentum
        dx = b * tpfl(10.0)  # initial separation
        
        # Initial positions (similar to setup_scattering)
        q1 = [-dx/ν2, -b/ν2, ν0]
        q2 = [dx/ν2, b/ν2, ν0]
        p1 = [p, ν0, ν0]
        p2 = [-p, ν0, ν0]
        
        # Calculate yba for both particles
        yba_12 = pomin.HamPM.ybaf(m[2], q2, q1, p2)  # particle 2 relative to 1
        yba_21 = pomin.HamPM.ybaf(m[1], q1, q2, p1)  # particle 1 relative to 2
        
        println("  yba_12 = $(yba_12)")
        println("  yba_21 = $(yba_21)")
        
        # Check for dangerous cancellation: yba + 1 ≈ 0
        cancellation_12 = yba_12 + ν1
        cancellation_21 = yba_21 + ν1
        
        println("  yba_12 + 1 = $(cancellation_12)")
        println("  yba_21 + 1 = $(cancellation_21)")
        
        # Check magnitude of cancellation terms
        if abs(cancellation_12) < tpfl(1e-10)
            println("  ⚠️  CRITICAL: yba_12 + 1 ≈ 0 (cancellation!)")
        end
        
        if abs(cancellation_21) < tpfl(1e-10)
            println("  ⚠️  CRITICAL: yba_21 + 1 ≈ 0 (cancellation!)")
        end
        
        # Test the H3 denominator term: (yba+1)^2
        denom_12 = (yba_12 + ν1)^2
        denom_21 = (yba_21 + ν1)^2
        
        println("  (yba_12 + 1)^2 = $(denom_12)")
        println("  (yba_21 + 1)^2 = $(denom_21)")
        
        # Check for very small denominators that could cause precision issues
        if denom_12 < tpfl(1e-20)
            println("  ⚠️  CRITICAL: (yba_12 + 1)^2 very small - potential precision explosion!")
        end
        
        if denom_21 < tpfl(1e-20)
            println("  ⚠️  CRITICAL: (yba_21 + 1)^2 very small - potential precision explosion!")
        end
        
        # Test the full H3 coefficient
        rab = pomin.HamPM.rf(q1, q2)
        Ena = pomin.HamPM.Enf(m[1], pomin.HamPM.psf(p1))
        Enb = pomin.HamPM.Enf(m[2], pomin.HamPM.psf(p2))
        
        if abs(denom_12) > tpfl(1e-30)
            h3_coeff_12 = ν1/(ν4*rab*Ena*Enb*yba_12*denom_12)
            println("  H3 coefficient (1->2): $(h3_coeff_12)")
            
            if abs(h3_coeff_12) > tpfl(1e10)
                println("  ⚠️  WARNING: H3 coefficient very large - potential precision issue!")
            end
        end
        
        if abs(denom_21) > tpfl(1e-30)
            h3_coeff_21 = ν1/(ν4*rab*Enb*Ena*yba_21*denom_21)
            println("  H3 coefficient (2->1): $(h3_coeff_21)")
            
            if abs(h3_coeff_21) > tpfl(1e10)
                println("  ⚠️  WARNING: H3 coefficient very large - potential precision issue!")
            end
        end
    end
end

function test_sqrt_precision_in_yba()
    println("\n=== Testing sqrt Precision in yba Calculation ===")
    
    setprecision(BigFloat, 128)
    tpfl = BigFloat
    
    # Test case where sqrt argument might be close to zero
    mb = tpfl(1.0)
    
    # Test different Θab values
    theta_values = [tpfl(0.0), tpfl(1e-10), tpfl(1e-5), tpfl(1e-2), tpfl(1.0)]
    
    for theta in theta_values
        println("\nTesting with Θab = $(theta)")
        
        # Calculate components of yba
        sqrt_arg = mb^2 + theta^2
        sqrt_result = sqrt(sqrt_arg)
        
        println("  mb^2 + Θab^2 = $(sqrt_arg)")
        println("  sqrt(mb^2 + Θab^2) = $(sqrt_result)")
        
        # Check for precision loss in sqrt
        expected_sqrt = mb * sqrt(tpfl(1) + (theta/mb)^2)  # Mathematically equivalent
        sqrt_diff = abs(sqrt_result - expected_sqrt)
        
        println("  Alternative calculation: $(expected_sqrt)")
        println("  Difference: $(sqrt_diff)")
        
        if sqrt_diff > tpfl(1e-25)
            println("  ⚠️  Potential sqrt precision issue!")
        else
            println("  ✅ sqrt precision OK")
        end
    end
end

function main()
    println("Square Root Cancellation Diagnostics")
    println("=" ^ 50)
    
    # Test 1: yba cancellation issues
    test_yba_cancellation()
    
    # Test 2: sqrt precision in yba calculation
    test_sqrt_precision_in_yba()
    
    println("\n" * ("=" ^ 50))
    println("Square Root Cancellation Diagnostics Complete")
    
    println("\n=== SUMMARY ===")
    println("This script tests for numerical cancellation in:")
    println("1. yba + 1 ≈ 0 (catastrophic cancellation)")
    println("2. (yba + 1)^2 in H3 denominators")
    println("3. sqrt precision in yba calculations")
    println("\nLook for WARNING messages indicating potential precision issues.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
