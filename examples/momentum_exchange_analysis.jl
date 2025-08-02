#!/usr/bin/env julia

"""
Momentum Exchange Analysis Script

This script analyzes momentum exchange in two-particle scattering across different
impact parameters and generates data suitable for log-log plotting.

Based on the MXIC.jl pattern but focused on systematic analysis.
"""

using LinearAlgebra
using DoubleFloats
using CSV

# Include necessary files from the project
include("../src/core/pomin-types.jl")
include("../src/physics/initial_data/idgen.jl")

"""
    analytical_momentum_change(p, b, d, m1, m2; G=1, c=1)

Calculate the analytical change in momentum for a two-particle scattering event.
Based on post-Minkowskian theory for gravitational scattering.
"""
function analytical_momentum_change(p, b, d, m1, m2; G=1, c=1)
    # Energy calculations for relativistic particles
    E1 = c * sqrt((c*m1)^2 + p^2)
    E2 = c * sqrt((c*m2)^2 + p^2)
    
    # Post-Minkowskian momentum change formula
    # This is a simplified version - the full formula would include higher-order terms
    dp = (2*G*(E1*E2)^2)/(b*p*(E1+E2)) * 
         (1 + (1/E1^2 + 1/E2^2 + 4/(E1*E2))*p^2 + p^4/(E1*E2)^2)
    
    return dp
end

"""
    analyze_momentum_exchange_vs_impact_parameter(m1, m2, p, d; 
                                                  b_min=0.01, b_max=10.0, n_points=50,
                                                  tpfl=Float64, save_csv=true)

Analyze momentum exchange across a range of impact parameters.
Returns data suitable for log-log plotting.
"""
function analyze_momentum_exchange_vs_impact_parameter(m1, m2, p, d; 
                                                       b_min=0.01, b_max=10.0, n_points=50,
                                                       tpfl=Float64, save_csv=true)
    
    # Convert parameters to specified type
    m1, m2, p, d = tpfl(m1), tpfl(m2), tpfl(p), tpfl(d)
    b_min, b_max = tpfl(b_min), tpfl(b_max)
    
    # Generate logarithmically spaced impact parameters
    log_b_min = log10(Float64(b_min))
    log_b_max = log10(Float64(b_max))
    b_range = [tpfl(10^x) for x in range(log_b_min, log_b_max, length=n_points)]
    
    # Storage for results
    analytical_results = tpfl[]
    
    println("Analyzing momentum exchange for $n_points impact parameters...")
    println("Mass 1: $m1, Mass 2: $m2")
    println("Momentum: $p, Initial separation: $d")
    println("Impact parameter range: $b_min to $b_max")
    println()
    
    # Calculate analytical results for each impact parameter
    for (i, b) in enumerate(b_range)
        if i % 10 == 0 || i == 1
            print("Progress: $i/$n_points (b = $(Float64(b)))\r")
        end
        
        # Calculate analytical momentum change
        dp_analytical = analytical_momentum_change(p, b, d, m1, m2)
        push!(analytical_results, dp_analytical)
    end
    println("\nCalculation complete!")
    
    # Create results dictionary
    results = (
        impact_parameters = b_range,
        analytical_momentum_change = analytical_results,
        parameters = (m1=m1, m2=m2, p=p, d=d),
        scaling_info = analyze_scaling_behavior(b_range, analytical_results)
    )
    
    # Save to CSV if requested
    if save_csv
        filename = "momentum_exchange_analysis_m1_$(Float64(m1))_m2_$(Float64(m2))_p_$(Float64(p)).csv"
        csv_data = Dict(
            "impact_parameter" => Float64.(b_range),
            "analytical_momentum_change" => Float64.(analytical_results),
            "log10_impact_parameter" => log10.(Float64.(b_range)),
            "log10_momentum_change" => log10.(Float64.(analytical_results))
        )
        
        CSV.write(filename, csv_data)
        println("Results saved to: $filename")
    end
    
    return results
end

"""
    analyze_scaling_behavior(b_range, dp_values)

Analyze the scaling behavior of momentum change vs impact parameter.
"""
function analyze_scaling_behavior(b_range, dp_values)
    n = length(b_range)
    if n < 3
        return (slope=NaN, correlation=NaN, scaling_regime="insufficient_data")
    end
    
    # Convert to log scale
    log_b = log10.(Float64.(b_range))
    log_dp = log10.(Float64.(dp_values))
    
    # Simple linear regression in log-log space
    n_fit = Float64(n)
    sum_log_b = sum(log_b)
    sum_log_dp = sum(log_dp)
    sum_log_b_sq = sum(log_b.^2)
    sum_log_b_log_dp = sum(log_b .* log_dp)
    
    # Calculate slope (scaling exponent)
    slope = (n_fit * sum_log_b_log_dp - sum_log_b * sum_log_dp) / 
            (n_fit * sum_log_b_sq - sum_log_b^2)
    
    # Calculate correlation coefficient
    mean_log_b = sum_log_b / n_fit
    mean_log_dp = sum_log_dp / n_fit
    
    ss_tot = sum((log_dp .- mean_log_dp).^2)
    ss_res = sum((log_dp .- (slope .* log_b .+ (mean_log_dp - slope * mean_log_b))).^2)
    r_squared = 1 - ss_res / ss_tot
    
    # Determine scaling regime
    scaling_regime = if abs(slope + 1.0) < 0.1
        "inverse_proportional"  # dp ∝ 1/b
    elseif abs(slope + 2.0) < 0.1
        "inverse_square"        # dp ∝ 1/b²
    elseif abs(slope) < 0.1
        "constant"              # dp ≈ constant
    else
        "power_law_$(round(slope, digits=2))"
    end
    
    return (slope=slope, r_squared=r_squared, scaling_regime=scaling_regime)
end

"""
    print_analysis_summary(results)

Print a summary of the momentum exchange analysis.
"""
function print_analysis_summary(results)
    println("\n" * "="^60)
    println("MOMENTUM EXCHANGE ANALYSIS SUMMARY")
    println("="^60)
    
    params = results.parameters
    println("Physical Parameters:")
    println("  Mass 1: $(params.m1)")
    println("  Mass 2: $(params.m2)")
    println("  Momentum: $(params.p)")
    println("  Initial separation: $(params.d)")
    
    b_range = results.impact_parameters
    dp_range = results.analytical_momentum_change
    
    println("\nImpact Parameter Range:")
    println("  Minimum: $(Float64(minimum(b_range)))")
    println("  Maximum: $(Float64(maximum(b_range)))")
    println("  Number of points: $(length(b_range))")
    
    println("\nMomentum Change Range:")
    println("  Minimum: $(Float64(minimum(dp_range)))")
    println("  Maximum: $(Float64(maximum(dp_range)))")
    println("  Dynamic range: $(Float64(maximum(dp_range)/minimum(dp_range)))")
    
    scaling = results.scaling_info
    println("\nScaling Analysis:")
    println("  Power law exponent: $(round(scaling.slope, digits=3))")
    println("  R² correlation: $(round(scaling.r_squared, digits=4))")
    println("  Scaling regime: $(scaling.scaling_regime)")
    
    println("\nFor log-log plotting:")
    println("  X-axis: log₁₀(impact parameter)")
    println("  Y-axis: log₁₀(momentum change)")
    println("  Expected slope: $(round(scaling.slope, digits=2))")
    
    println("="^60)
end

"""
    run_comprehensive_analysis()

Run a comprehensive momentum exchange analysis with multiple parameter sets.
"""
function run_comprehensive_analysis()
    println("COMPREHENSIVE MOMENTUM EXCHANGE ANALYSIS")
    println("="^50)
    
    # Test cases
    test_cases = [
        (m1=0.1, m2=0.1, p=1.0, d=10.0, name="equal_mass_moderate"),
        (m1=0.01, m2=0.01, p=0.5, d=20.0, name="equal_mass_low_momentum"),
        (m1=1.0, m2=0.1, p=1.0, d=10.0, name="mass_ratio_10_to_1"),
        (m1=0.0, m2=0.0, p=1.0, d=10.0, name="massless_particles")
    ]
    
    all_results = []
    
    for (i, case) in enumerate(test_cases)
        println("\n" * "-"^50)
        println("Test Case $i: $(case.name)")
        println("-"^50)
        
        try
            results = analyze_momentum_exchange_vs_impact_parameter(
                case.m1, case.m2, case.p, case.d;
                b_min=0.01, b_max=10.0, n_points=30,
                save_csv=true
            )
            
            print_analysis_summary(results)
            push!(all_results, (case=case, results=results))
            
        catch e
            println("Error in test case $(case.name): $e")
        end
    end
    
    println("\n" * "="^50)
    println("ANALYSIS COMPLETE")
    println("Generated CSV files can be plotted with external tools.")
    println("For log-log plots, use log10(impact_parameter) vs log10(momentum_change)")
    println("="^50)
    
    return all_results
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("Starting momentum exchange analysis...")
    
    # Run single analysis example
    println("\n1. Single Analysis Example:")
    results = analyze_momentum_exchange_vs_impact_parameter(
        0.1, 0.1, 1.0, 10.0;  # m1, m2, p, d
        b_min=0.01, b_max=10.0, n_points=25
    )
    print_analysis_summary(results)
    
    # Run comprehensive analysis
    println("\n\n2. Comprehensive Analysis:")
    all_results = run_comprehensive_analysis()
    
    println("\nAnalysis complete! Check the generated CSV files for plotting data.")
end
