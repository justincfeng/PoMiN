#!/usr/bin/env julia

"""
Simple Momentum Exchange Analysis

This script analyzes momentum exchange in two-particle scattering across different
impact parameters, similar to MXIC.jl but focused on generating log plot data.

Uses only basic Julia functionality - no external dependencies.
"""

using LinearAlgebra

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
    generate_log_plot_data(m1, m2, p, d; b_min=0.01, b_max=10.0, n_points=50)

Generate data for log-log plot of momentum exchange vs impact parameter.
Returns arrays suitable for plotting.
"""
function generate_log_plot_data(m1, m2, p, d; b_min=0.01, b_max=10.0, n_points=50)
    
    # Generate logarithmically spaced impact parameters
    log_b_min = log10(b_min)
    log_b_max = log10(b_max)
    b_range = [10^x for x in range(log_b_min, log_b_max, length=n_points)]
    
    # Calculate analytical momentum change for each impact parameter
    dp_analytical = Float64[]
    
    println("Generating momentum exchange data...")
    println("Parameters: m1=$m1, m2=$m2, p=$p, d=$d")
    println("Impact parameter range: $b_min to $b_max ($n_points points)")
    
    for (i, b) in enumerate(b_range)
        dp = analytical_momentum_change(p, b, d, m1, m2)
        push!(dp_analytical, dp)
        
        if i % 10 == 0 || i == 1
            println("  Point $i/$n_points: b=$b, Δp=$dp")
        end
    end
    
    return b_range, dp_analytical
end

"""
    save_data_for_plotting(filename, b_range, dp_values; m1, m2, p, d)

Save data in a simple format for external plotting.
"""
function save_data_for_plotting(filename, b_range, dp_values; m1, m2, p, d)
    open(filename, "w") do file
        # Write header with parameters
        println(file, "# Momentum Exchange Analysis Data")
        println(file, "# Parameters: m1=$m1, m2=$m2, p=$p, d=$d")
        println(file, "# Columns: impact_parameter, momentum_change, log10_b, log10_dp")
        
        # Write data
        for (b, dp) in zip(b_range, dp_values)
            log_b = log10(b)
            log_dp = log10(dp)
            println(file, "$b,$dp,$log_b,$log_dp")
        end
    end
    println("Data saved to: $filename")
end

"""
    analyze_scaling(b_range, dp_values)

Analyze the scaling behavior in log-log space.
"""
function analyze_scaling(b_range, dp_values)
    log_b = log10.(b_range)
    log_dp = log10.(dp_values)
    
    # Simple linear regression in log-log space
    n = length(log_b)
    sum_log_b = sum(log_b)
    sum_log_dp = sum(log_dp)
    sum_log_b_sq = sum(log_b.^2)
    sum_log_b_log_dp = sum(log_b .* log_dp)
    
    # Calculate slope (scaling exponent)
    slope = (n * sum_log_b_log_dp - sum_log_b * sum_log_dp) / 
            (n * sum_log_b_sq - sum_log_b^2)
    
    # Calculate R²
    mean_log_b = sum_log_b / n
    mean_log_dp = sum_log_dp / n
    
    ss_tot = sum((log_dp .- mean_log_dp).^2)
    predicted = slope .* log_b .+ (mean_log_dp - slope * mean_log_b)
    ss_res = sum((log_dp .- predicted).^2)
    r_squared = 1 - ss_res / ss_tot
    
    return slope, r_squared
end

"""
    print_analysis_summary(b_range, dp_values, m1, m2, p, d)

Print a summary of the analysis results.
"""
function print_analysis_summary(b_range, dp_values, m1, m2, p, d)
    slope, r_squared = analyze_scaling(b_range, dp_values)
    
    println("\n" * "="^60)
    println("MOMENTUM EXCHANGE ANALYSIS SUMMARY")
    println("="^60)
    
    println("Physical Parameters:")
    println("  Mass 1: $m1")
    println("  Mass 2: $m2") 
    println("  Momentum: $p")
    println("  Initial separation: $d")
    
    println("\nData Range:")
    println("  Impact parameters: $(minimum(b_range)) to $(maximum(b_range))")
    println("  Momentum changes: $(minimum(dp_values)) to $(maximum(dp_values))")
    println("  Dynamic range: $(maximum(dp_values)/minimum(dp_values))")
    
    println("\nScaling Analysis (log-log):")
    println("  Power law exponent: $(round(slope, digits=3))")
    println("  R² correlation: $(round(r_squared, digits=4))")
    
    expected_scaling = if abs(slope + 1.0) < 0.1
        "≈ 1/b (inverse proportional)"
    elseif abs(slope + 2.0) < 0.1
        "≈ 1/b² (inverse square)"
    else
        "power law with exponent $(round(slope, digits=2))"
    end
    println("  Scaling behavior: $expected_scaling")
    
    println("\nFor Log-Log Plotting:")
    println("  X-axis: log₁₀(impact parameter)")
    println("  Y-axis: log₁₀(momentum change)")
    println("  Expected slope: $(round(slope, digits=2))")
    
    println("="^60)
end

# Main analysis function
function run_momentum_exchange_analysis()
    println("MOMENTUM EXCHANGE ANALYSIS")
    println("="^40)
    
    # Test case 1: Equal mass particles
    println("\n1. Equal Mass Particles (m1 = m2 = 0.1)")
    b_range1, dp_values1 = generate_log_plot_data(0.1, 0.1, 1.0, 10.0; n_points=30)
    save_data_for_plotting("momentum_exchange_equal_mass.dat", b_range1, dp_values1; 
                          m1=0.1, m2=0.1, p=1.0, d=10.0)
    print_analysis_summary(b_range1, dp_values1, 0.1, 0.1, 1.0, 10.0)
    
    # Test case 2: Mass ratio 10:1
    println("\n2. Unequal Masses (m1 = 1.0, m2 = 0.1)")
    b_range2, dp_values2 = generate_log_plot_data(1.0, 0.1, 1.0, 10.0; n_points=30)
    save_data_for_plotting("momentum_exchange_mass_ratio.dat", b_range2, dp_values2;
                          m1=1.0, m2=0.1, p=1.0, d=10.0)
    print_analysis_summary(b_range2, dp_values2, 1.0, 0.1, 1.0, 10.0)
    
    # Test case 3: Massless particles
    println("\n3. Massless Particles (m1 = m2 = 0.0)")
    b_range3, dp_values3 = generate_log_plot_data(0.0, 0.0, 1.0, 10.0; n_points=30)
    save_data_for_plotting("momentum_exchange_massless.dat", b_range3, dp_values3;
                          m1=0.0, m2=0.0, p=1.0, d=10.0)
    print_analysis_summary(b_range3, dp_values3, 0.0, 0.0, 1.0, 10.0)
    
    println("\n" * "="^40)
    println("ANALYSIS COMPLETE")
    println("Generated data files:")
    println("- momentum_exchange_equal_mass.dat")
    println("- momentum_exchange_mass_ratio.dat") 
    println("- momentum_exchange_massless.dat")
    println("\nTo create log plots:")
    println("- Plot column 3 (log10_b) vs column 4 (log10_dp)")
    println("- Or use any plotting tool with columns 1 and 2 on log scales")
    println("="^40)
end

# Run the analysis if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_momentum_exchange_analysis()
end
