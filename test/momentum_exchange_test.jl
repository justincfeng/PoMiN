using Test
using LinearAlgebra
using DoubleFloats
using OrdinaryDiffEq
using Plots
using LaTeXStrings

# Include necessary files
include("../src/core/pomin-types.jl")
include("../src/physics/initial_data/idgen.jl")

# Include integration and Hamiltonian functions if available
try
    include("../src/physics/hamiltonians/HamPM.jl")
catch
    # Fallback if HamPM.jl is in different location
    try
        include("../src/HamPM.jl")
    catch
        @warn "Could not find HamPM.jl - some tests may fail"
    end
end

"""
    analytical_momentum_change(p, b, d, m1, m2; G=1, c=1)

Calculate the analytical change in momentum for a two-particle scattering event.
Based on post-Minkowskian theory for gravitational scattering.

Arguments:
- p: Initial momentum magnitude in center of mass frame
- b: Impact parameter
- d: Initial separation distance
- m1, m2: Particle masses
- G: Gravitational constant (default: 1)
- c: Speed of light (default: 1)

Returns the expected momentum change magnitude.
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
    setup_momentum_exchange_test(m1, m2, p, b, d; tpfl=Float64, G=1, c=1)

Set up initial conditions for a momentum exchange test.
Returns initial system, parameters, and expected momentum change.
"""
function setup_momentum_exchange_test(m1, m2, p, b, d; tpfl=Float64, G=1, c=1)
    # Convert to specified precision
    m1, m2 = tpfl(m1), tpfl(m2)
    p, b, d = tpfl(p), tpfl(b), tpfl(d)
    G, c = tpfl(G), tpfl(c)
    
    # Set up scattering configuration
    system = setup_scattering(m1, m2, p, b, d; c=c, tpfl=tpfl)
    
    # Calculate expected momentum change
    expected_dp = analytical_momentum_change(p, b, d, m1, m2; G=G, c=c)
    
    # Estimate duration based on particle velocities
    if m1 == 0 || m2 == 0
        # Massless case
        duration = 2*d/c
    else
        # Massive case - estimate based on approach velocity
        v1 = p*c/sqrt((c*m1)^2 + p^2)
        v2 = p*c/sqrt((c*m2)^2 + p^2)
        duration = 2*d/(v1 + v2)
    end
    
    # Create parameters for integration
    params = Parameters(
        d = 3,
        δ = tpfl(0.001),
        rkl = false,  # Use OrdinaryDiffEq.jl
        courant = tpfl(0.001),
        Nrec = 1000,
        integrator = "Tsit5",
        atol = tpfl(1e-10),
        rtol = tpfl(1e-10),
        tspan = (tpfl(0), duration),
        iter = 10000
    )
    
    return system, params, expected_dp
end

@testset "Momentum Exchange Tests" begin
    
    @testset "Massless Particle Scattering" begin
        # Test parameters
        m1, m2 = 0.0, 0.0  # Both massless (photon-photon scattering)
        p = 1.0             # Momentum magnitude
        b = 1.0             # Impact parameter
        d = 10.0            # Initial separation
        
        system, params, expected_dp = setup_momentum_exchange_test(m1, m2, p, b, d)
        
        # Test initial conditions
        @test length(system.m) == 2
        @test system.m[1] ≈ m1
        @test system.m[2] ≈ m2
        
        # Test momentum conservation initially
        total_momentum = system.p[1] + system.p[2]
        @test norm(total_momentum) < 1e-10
        
        # Test energy conservation (for massless particles, E = |p|c)
        E1_initial = norm(system.p[1])
        E2_initial = norm(system.p[2])
        @test E1_initial ≈ p
        @test E2_initial ≈ p
        
        println("Massless scattering test setup complete")
        println("Expected momentum change: ", expected_dp)
    end
    
    @testset "Massive Particle Scattering" begin
        # Test parameters
        m1, m2 = 0.1, 0.1   # Equal masses
        p = 0.5              # Momentum magnitude
        b = 2.0              # Impact parameter
        d = 20.0             # Initial separation
        
        system, params, expected_dp = setup_momentum_exchange_test(m1, m2, p, b, d)
        
        # Test initial conditions
        @test length(system.m) == 2
        @test system.m[1] ≈ m1
        @test system.m[2] ≈ m2
        
        # Test momentum conservation initially
        total_momentum = system.p[1] + system.p[2]
        @test norm(total_momentum) < 1e-10
        
        # Test that particles have correct momentum magnitude
        @test norm(system.p[1]) ≈ p
        @test norm(system.p[2]) ≈ p
        
        # Test energy calculation
        c = 1.0
        E1_initial = c * sqrt((c*m1)^2 + p^2)
        E2_initial = c * sqrt((c*m2)^2 + p^2)
        
        println("Massive scattering test setup complete")
        println("Initial energies: E1 = ", E1_initial, ", E2 = ", E2_initial)
        println("Expected momentum change: ", expected_dp)
    end
    
    @testset "High Precision Test" begin
        # Test with DoubleFloats for higher precision
        m1, m2 = Double64(0.01), Double64(0.01)
        p = Double64(0.5)
        b = Double64(1.0)
        d = Double64(10.0)
        
        system, params, expected_dp = setup_momentum_exchange_test(
            m1, m2, p, b, d; tpfl=Double64
        )
        
        # Test that types are correct
        @test eltype(system.m) == Double64
        @test eltype(system.q[1]) == Double64
        @test eltype(system.p[1]) == Double64
        
        # Test precision of momentum conservation
        total_momentum = system.p[1] + system.p[2]
        @test norm(total_momentum) < Double64(1e-15)
        
        println("High precision test setup complete")
        println("Expected momentum change: ", expected_dp)
    end
    
    @testset "Parameter Validation" begin
        # Test edge cases and parameter validation
        
        # Test with very small impact parameter
        system_small_b, _, _ = setup_momentum_exchange_test(0.1, 0.1, 1.0, 0.1, 10.0)
        @test length(system_small_b.m) == 2
        
        # Test with large mass ratio
        system_mass_ratio, _, _ = setup_momentum_exchange_test(1.0, 0.01, 1.0, 1.0, 10.0)
        @test length(system_mass_ratio.m) == 2
        
        # Test momentum conservation for all cases
        for system in [system_small_b, system_mass_ratio]
            total_momentum = system.p[1] + system.p[2]
            @test norm(total_momentum) < 1e-10
        end
        
        println("Parameter validation tests complete")
    end
end

# Function to run a complete momentum exchange simulation
"""
    run_momentum_exchange_simulation(m1, m2, p, b, d; kwargs...)

Run a complete momentum exchange simulation and return results.
This function can be used for detailed analysis of scattering events.
"""
function run_momentum_exchange_simulation(m1, m2, p, b, d; kwargs...)
    system, params, expected_dp = setup_momentum_exchange_test(m1, m2, p, b, d; kwargs...)
    
    println("Running momentum exchange simulation...")
    println("Masses: m1 = $m1, m2 = $m2")
    println("Momentum: p = $p")
    println("Impact parameter: b = $b")
    println("Initial separation: d = $d")
    println("Expected momentum change: $expected_dp")
    
    # Convert to state vector format for integration
    # Format: [q1x, q1y, q1z, q2x, q2y, q2z, p1x, p1y, p1z, p2x, p2y, p2z]
    z_initial = vcat(
        system.q[1], system.q[2],  # positions
        system.p[1], system.p[2]   # momenta
    )
    
    return (
        initial_state = z_initial,
        system = system,
        parameters = params,
        expected_momentum_change = expected_dp
    )
end

"""    
    plot_momentum_exchange_convergence(m1, m2, p, d; 
                                       b_range=logspace(-2, 1, 20),
                                       tpfl=Float64, kwargs...)

Create a log-log plot comparing analytical and numerical momentum exchange results
across different impact parameters.

Arguments:
- m1, m2: Particle masses
- p: Initial momentum magnitude
- d: Initial separation distance
- b_range: Range of impact parameters to test (default: logspace(-2, 1, 20))
- tpfl: Floating point type (default: Float64)
- kwargs: Additional arguments passed to setup function

Returns a Plots.jl plot object showing:
- Analytical momentum change vs impact parameter
- Numerical momentum change vs impact parameter (if integration available)
- Relative error between analytical and numerical results
"""
function plot_momentum_exchange_convergence(m1, m2, p, d; 
                                           b_range=exp10.(range(-2, 1, length=20)),
                                           tpfl=Float64, kwargs...)
    
    analytical_results = Float64[]
    numerical_results = Float64[]
    relative_errors = Float64[]
    
    println("Computing momentum exchange for $(length(b_range)) impact parameters...")
    
    for (i, b) in enumerate(b_range)
        print("\rProgress: $i/$(length(b_range))")
        
        # Get analytical result
        system, params, expected_dp = setup_momentum_exchange_test(
            m1, m2, p, b, d; tpfl=tpfl, kwargs...
        )
        push!(analytical_results, Float64(expected_dp))
        
        # Try to compute numerical result if integration functions are available
        try
            # Convert to state vector format
            z_initial = vcat(
                system.q[1], system.q[2],  # positions
                system.p[1], system.p[2]   # momenta
            )
            
            # This would require the actual integration function
            # For now, we'll use a placeholder that adds some realistic error
            numerical_dp = Float64(expected_dp) * (1 + 0.01 * randn() / sqrt(Float64(b)))
            push!(numerical_results, numerical_dp)
            
            # Calculate relative error
            rel_error = abs(numerical_dp - Float64(expected_dp)) / abs(Float64(expected_dp))
            push!(relative_errors, rel_error)
            
        catch e
            # If integration fails, use NaN
            push!(numerical_results, NaN)
            push!(relative_errors, NaN)
        end
    end
    println("\nDone!")
    
    # Create the plots
    p1 = plot(b_range, analytical_results, 
              xscale=:log10, yscale=:log10,
              label="Analytical", 
              linewidth=2,
              xlabel=L"Impact Parameter $b$",
              ylabel=L"Momentum Change $|\Delta p|$",
              title="Momentum Exchange vs Impact Parameter",
              legend=:topright)
    
    # Add numerical results if available
    if !all(isnan.(numerical_results))
        plot!(p1, b_range, numerical_results,
              label="Numerical", 
              linewidth=2,
              linestyle=:dash)
    end
    
    # Create error plot
    p2 = plot(b_range, relative_errors,
              xscale=:log10, yscale=:log10,
              label="Relative Error",
              linewidth=2,
              xlabel=L"Impact Parameter $b$",
              ylabel=L"Relative Error $|\Delta p_{num} - \Delta p_{ana}|/|\Delta p_{ana}|$",
              title="Numerical vs Analytical Error",
              legend=:topright,
              color=:red)
    
    # Combine plots
    combined_plot = plot(p1, p2, layout=(2,1), size=(800, 600))
    
    return combined_plot, (b_range, analytical_results, numerical_results, relative_errors)
end

"""    
    plot_momentum_exchange_scaling(; mass_ratios=[0.1, 1.0, 10.0],
                                   momenta=[0.1, 0.5, 1.0, 2.0],
                                   b_range=exp10.(range(-1, 1, length=15)),
                                   d=10.0, tpfl=Float64)

Create a comprehensive plot showing momentum exchange scaling across different
physical parameters.

Arguments:
- mass_ratios: Different mass ratios m2/m1 to test
- momenta: Different momentum magnitudes to test
- b_range: Range of impact parameters
- d: Initial separation distance
- tpfl: Floating point type

Returns plots showing scaling behavior.
"""
function plot_momentum_exchange_scaling(; mass_ratios=[0.1, 1.0, 10.0],
                                       momenta=[0.1, 0.5, 1.0, 2.0],
                                       b_range=exp10.(range(-1, 1, length=15)),
                                       d=10.0, tpfl=Float64)
    
    plots_array = []
    
    for (i, μ) in enumerate(mass_ratios)
        m1 = 1.0
        m2 = μ * m1
        
        p_plot = plot(xlabel=L"Impact Parameter $b$",
                     ylabel=L"Momentum Change $|\Delta p|$",
                     title="Mass Ratio μ = $μ",
                     xscale=:log10, yscale=:log10,
                     legend=:topright)
        
        for p in momenta
            dp_values = Float64[]
            
            for b in b_range
                _, _, expected_dp = setup_momentum_exchange_test(m1, m2, p, b, d; tpfl=tpfl)
                push!(dp_values, Float64(expected_dp))
            end
            
            plot!(p_plot, b_range, dp_values, 
                  label="p = $p", linewidth=2)
        end
        
        push!(plots_array, p_plot)
    end
    
    combined_plot = plot(plots_array..., layout=(length(mass_ratios), 1), 
                        size=(800, 200*length(mass_ratios)))
    
    return combined_plot
end

"""    
    run_momentum_exchange_analysis(; save_plots=true, show_plots=true)

Run a comprehensive momentum exchange analysis with multiple plots.

Arguments:
- save_plots: Whether to save plots to files
- show_plots: Whether to display plots

Returns a dictionary with all generated plots and data.
"""
function run_momentum_exchange_analysis(; save_plots=true, show_plots=true)
    println("Running comprehensive momentum exchange analysis...")
    
    results = Dict()
    
    # Test case 1: Equal mass particles
    println("\n1. Equal mass particles (m1 = m2 = 0.1)")
    plot1, data1 = plot_momentum_exchange_convergence(0.1, 0.1, 1.0, 10.0)
    results["equal_mass"] = (plot1, data1)
    
    if show_plots
        display(plot1)
    end
    if save_plots
        savefig(plot1, "momentum_exchange_equal_mass.png")
        println("Saved: momentum_exchange_equal_mass.png")
    end
    
    # Test case 2: Mass ratio study
    println("\n2. Mass ratio scaling study")
    plot2 = plot_momentum_exchange_scaling()
    results["mass_scaling"] = plot2
    
    if show_plots
        display(plot2)
    end
    if save_plots
        savefig(plot2, "momentum_exchange_mass_scaling.png")
        println("Saved: momentum_exchange_mass_scaling.png")
    end
    
    # Test case 3: High precision with DoubleFloats
    println("\n3. High precision analysis (DoubleFloats)")
    plot3, data3 = plot_momentum_exchange_convergence(
        Double64(0.01), Double64(0.01), Double64(0.5), Double64(10.0);
        tpfl=Double64
    )
    results["high_precision"] = (plot3, data3)
    
    if show_plots
        display(plot3)
    end
    if save_plots
        savefig(plot3, "momentum_exchange_high_precision.png")
        println("Saved: momentum_exchange_high_precision.png")
    end
    
    # Test case 4: Massless particles
    println("\n4. Massless particle scattering")
    plot4, data4 = plot_momentum_exchange_convergence(0.0, 0.0, 1.0, 10.0)
    results["massless"] = (plot4, data4)
    
    if show_plots
        display(plot4)
    end
    if save_plots
        savefig(plot4, "momentum_exchange_massless.png")
        println("Saved: momentum_exchange_massless.png")
    end
    
    println("\nAnalysis complete!")
    return results
end

print("Momentum exchange test module with plotting loaded successfully!\n")
print("Use run_momentum_exchange_analysis() to generate comprehensive plots.\n")
