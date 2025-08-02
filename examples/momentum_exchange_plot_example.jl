#!/usr/bin/env julia

"""
Example script demonstrating momentum exchange plotting functionality.

This script shows how to use the momentum exchange test module to create
log plots comparing analytical and numerical results across different
impact parameters.
"""

# Add the test directory to the path so we can include the test module
push!(LOAD_PATH, "../test")

using Plots
using LaTeXStrings

# Include the momentum exchange test module
include("../test/momentum_exchange_test.jl")

println("=== Momentum Exchange Plotting Example ===\n")

# Example 1: Basic momentum exchange plot for equal mass particles
println("1. Creating basic momentum exchange plot...")
try
    plot1, data1 = plot_momentum_exchange_convergence(
        0.1,    # m1
        0.1,    # m2  
        1.0,    # momentum magnitude
        10.0;   # initial separation
        b_range=exp10.(range(-1.5, 0.5, length=15))
    )
    
    display(plot1)
    savefig(plot1, "example_momentum_exchange_basic.png")
    println("✓ Saved: example_momentum_exchange_basic.png")
    
    # Print some data
    b_range, analytical, numerical, errors = data1
    println("Impact parameter range: $(minimum(b_range)) to $(maximum(b_range))")
    println("Analytical momentum change range: $(minimum(analytical)) to $(maximum(analytical))")
    
catch e
    println("Error in basic plot: $e")
end

println()

# Example 2: Mass ratio scaling study
println("2. Creating mass ratio scaling plot...")
try
    plot2 = plot_momentum_exchange_scaling(
        mass_ratios=[0.1, 1.0, 10.0],
        momenta=[0.5, 1.0, 2.0],
        b_range=exp10.(range(-1, 1, length=12))
    )
    
    display(plot2)
    savefig(plot2, "example_momentum_exchange_scaling.png")
    println("✓ Saved: example_momentum_exchange_scaling.png")
    
catch e
    println("Error in scaling plot: $e")
end

println()

# Example 3: High precision comparison
println("3. Creating high precision comparison...")
try
    using DoubleFloats
    
    plot3, data3 = plot_momentum_exchange_convergence(
        Double64(0.05),  # m1
        Double64(0.05),  # m2
        Double64(0.8),   # momentum
        Double64(15.0);  # separation
        tpfl=Double64,
        b_range=exp10.(range(-1.2, 0.8, length=12))
    )
    
    display(plot3)
    savefig(plot3, "example_momentum_exchange_precision.png")
    println("✓ Saved: example_momentum_exchange_precision.png")
    
catch e
    println("Error in precision plot: $e")
end

println()

# Example 4: Custom single plot with specific parameters
println("4. Creating custom single parameter plot...")
try
    # Parameters for a specific scattering scenario
    m1, m2 = 0.2, 0.1  # Unequal masses
    p = 1.5             # High momentum
    d = 8.0             # Moderate separation
    
    # Generate impact parameter range
    b_values = exp10.(range(-1.8, 0.2, length=18))
    
    # Calculate analytical results
    analytical_dp = Float64[]
    for b in b_values
        _, _, dp = setup_momentum_exchange_test(m1, m2, p, b, d)
        push!(analytical_dp, Float64(dp))
    end
    
    # Create custom plot
    custom_plot = plot(
        b_values, analytical_dp,
        xscale=:log10, yscale=:log10,
        xlabel=L"Impact Parameter $b$",
        ylabel=L"Momentum Change $|\Delta p|$",
        title="Custom Momentum Exchange: m₁=$m1, m₂=$m2, p=$p",
        linewidth=3,
        color=:blue,
        label="Analytical Result",
        grid=true,
        gridwidth=1,
        gridcolor=:gray,
        gridalpha=0.3
    )
    
    # Add theoretical scaling line for comparison
    # For small impact parameters, expect dp ∝ 1/b behavior
    b_theory = b_values[b_values .< 0.5]
    dp_theory = analytical_dp[1] * (b_values[1] ./ b_theory)
    plot!(custom_plot, b_theory, dp_theory,
          linestyle=:dash, color=:red, linewidth=2,
          label=L"$\propto 1/b$ scaling")
    
    display(custom_plot)
    savefig(custom_plot, "example_momentum_exchange_custom.png")
    println("✓ Saved: example_momentum_exchange_custom.png")
    
    # Print scaling analysis
    println("Scaling analysis:")
    for i in 2:length(b_values)
        if i <= 5  # Only print first few points
            ratio_b = b_values[i] / b_values[i-1]
            ratio_dp = analytical_dp[i-1] / analytical_dp[i]  # Note: inverse because dp decreases as b increases
            println("  b ratio: $(round(ratio_b, digits=3)), dp ratio: $(round(ratio_dp, digits=3))")
        end
    end
    
catch e
    println("Error in custom plot: $e")
end

println("\n=== Example Complete ===")
println("Generated plots:")
println("- example_momentum_exchange_basic.png")
println("- example_momentum_exchange_scaling.png") 
println("- example_momentum_exchange_precision.png")
println("- example_momentum_exchange_custom.png")
println("\nTo run the full analysis suite, use: run_momentum_exchange_analysis()")
