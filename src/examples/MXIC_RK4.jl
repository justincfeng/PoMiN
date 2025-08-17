#!/usr/bin/env julia

"""
PoMiN Momentum Exchange Sweep using RK4 Integrator with Courant Number (tcour)
Based on the ApJ paper implementation with adaptive timestepping.
"""

using Plots
using Printf
include("../pomin.jl")

"""
    dpfunc_analytical(p, b, m1, m2; G=1, c=1)

Analytical momentum exchange formula for post-Minkowskian scattering.
Returns the change in momentum for particle 1 in the y-direction.
"""
function dpfunc_analytical(p, b, m1, m2; Ga=1, ca=1)
    tpfl = typeof(p)
    ν0, ν1, ν2, ν3, ν4, ν5, ν6, ν7, ν8, ν9 = pomin.tpnum(tpfl)
    c = tpfl(ca)
    G = tpfl(Ga)

    E1 = c*sqrt((c*m1)^2 + p^2)
    E2 = c*sqrt((c*m2)^2 + p^2)
    # Formula matches LSB_momX.txt exactly
    dp = (ν2*G*(E1*E2)^2)/(p*(E1+E2)) *
         (ν1 + (ν1/E1^2 + ν1/E2^2 + ν4/(E1*E2))*p^2 + p^4/(E1*E2)^2) / b
    return dp
end

function main()
    println("=== PoMiN Momentum Exchange Sweep (RK4 + Courant) ===")
    println("Using RK4 integrator with adaptive timestepping based on CFL condition")
    println("Parameters from ApJ paper implementation")
    
    # Set BigFloat precision to 128 bits (quad precision) to eliminate error plateau
    setprecision(BigFloat, 128)
    println("Precision: BigFloat with $(precision(BigFloat)) bits (quad precision)\n")
    
    # Physical parameters (using BigFloat with quad precision to eliminate error plateau)
    m1 = BigFloat(1.0)
    m2 = BigFloat(1.0)  
    p = BigFloat(10.0)
    
    # Impact parameter sweep parameters (all BigFloat for consistency)
    b_base = BigFloat(10.0)
    scale_factor = BigFloat("1.1108305558745590335689712446765042841434478759765625")  # From alternative MXIC
    nn_values = 50:50:350  # Range from 1 to 250 in steps of 50
    
    results = []
    
    for nn in nn_values
        println("Testing with nn = $nn")
        
        # Generate impact parameter with wide dynamic range (ensure BigFloat precision)
        b = b_base * scale_factor^nn
        dx = BigFloat(1e9)*b  # Optimal separation (not too large, not too small)
        
        println("  Impact parameter b = $(round(b, sigdigits=6))")
        println("  Initial separation dx = $(round(dx, sigdigits=6))")
        
        # Set up scattering system (using BigFloat type for quad precision)
        system = pomin.setup_scattering(m1, m2, p, b, dx, tpfl=BigFloat)
        
        # Calculate duration using relativistic velocity calculation
        v1 = p / sqrt(m1^2 + p^2)  # Relativistic velocity for particle 1
        v2 = p / sqrt(m2^2 + p^2)  # Relativistic velocity for particle 2
        τ = dx / (v1 + v2)         # Time for particles to meet
        t_flight = BigFloat(10.0) * τ       # Total scattering time (same as MXIC.jl)
        
        # RK4 integrator parameters from ApJ paper (all BigFloat for precision)
        # δ = initial timestep, courant = Courant number for CFL condition
        δ_initial = BigFloat(1.0)          # Initial timestep (reduced for higher accuracy)
        courant = BigFloat(0.01)          # Courant number (reduced further for maximum accuracy)
        
        # Solve with RK4 integrator using adaptive timestepping (BigFloat precision)
        sol = pomin.solve(system, pomin.ParametersRK4((BigFloat(0.0), t_flight), 
                         δ=δ_initial, courant=courant, Nrec=-1))
        
        # Extract momentum change (particle 1, y-component) - FIXED to match analytical formula
        # Phase space vector: [q1x, q1y, q1z, q2x, q2y, q2z, p1x, p1y, p1z, p2x, p2y, p2z]
        # The analytical formula gives y-component change for particle 1
        
        # Debug: Print phase space vector structure for first iteration
        if nn == nn_values[1]
            println("  DEBUG: Initial state vector length = $(length(sol.z[1]))")
            println("  DEBUG: Initial state = $(sol.z[1])")
            println("  DEBUG: Expected order: [q1x, q1y, q1z, q2x, q2y, q2z, p1x, p1y, p1z, p2x, p2y, p2z]")
        end
        
        # Extract positions and momenta for analysis
        q1_initial = [sol.z[1][1], sol.z[1][2], sol.z[1][3]]   # Initial position particle 1
        q2_initial = [sol.z[1][4], sol.z[1][5], sol.z[1][6]]   # Initial position particle 2
        q1_final = [sol.z[end][1], sol.z[end][2], sol.z[end][3]]   # Final position particle 1
        q2_final = [sol.z[end][4], sol.z[end][5], sol.z[end][6]]   # Final position particle 2
        
        # Calculate initial and final separations
        r_initial = sqrt(sum((q1_initial - q2_initial).^2))
        r_final = sqrt(sum((q1_final - q2_final).^2))
        
        p1y_initial = sol.z[1][8]   # p1y component (8th element)
        p1y_final = sol.z[end][8]   # p1y component at end
        dp_numerical = p1y_final - p1y_initial  # Y-component only, not magnitude
        
        # Calculate analytical prediction
        dp_analytical = dpfunc_analytical(p, b, m1, m2)
        
        # Calculate relative error (ensure BigFloat precision throughout)
        rel_error = abs(abs(dp_numerical) - abs(dp_analytical)) / abs(dp_analytical) * BigFloat(100)
        
        println("  Analytical Δp = $(round(dp_analytical, sigdigits=8))")
        println("  Numerical  Δp = $(round(dp_numerical, sigdigits=8))")
        println("  Relative error = $(round(rel_error, sigdigits=4))%")
        println("  Integration time = $(round(t_flight, sigdigits=6))")
        println("  Initial separation = $(round(r_initial, sigdigits=6))")
        println("  Final separation = $(round(r_final, sigdigits=6))")
        println("  Separation ratio (final/initial) = $(round(r_final/r_initial, sigdigits=4))")
        println("  Saved $(length(sol.z)) time points")
        println()
        
        # Store results
        push!(results, (nn, b, dp_analytical, dp_numerical, rel_error))
    end
    
    # Print summary
    println("=== Summary ===")
    println("nn      b               Δp_analytical   Δp_numerical    Rel_Error(%)")
    println("----------------------------------------------------------------------")
    for (nn, b, dp_a, dp_n, err) in results
        @printf("%-8d%-16g%-16g%-16g%-8.2f\n", nn, b, dp_a, dp_n, err)
    end
    
    return results
end

"""
    plot_momentum_exchange_rk4(results)

Generate plots showing momentum exchange scaling and accuracy for RK4 integrator.
Matches the style and quality of MXIC.jl plots.
"""
function plot_momentum_exchange_rk4(results)
    # Extract data for plotting (convert BigFloat to Float64 for plotting)
    b_values = [Float64(r[2]) for r in results]
    dp_analytical = [Float64(r[3]) for r in results]
    dp_numerical = [Float64(r[4]) for r in results]
    rel_errors = [Float64(r[5]) for r in results]
    
    # Create log-log plot of momentum exchange vs impact parameter
    p1 = plot(b_values, dp_analytical, 
              xscale=:log10, yscale=:log10,
              marker=:circle, markersize=6, linewidth=2,
              label="Analytical (PM Theory)",
              xlabel="Impact Parameter b", ylabel="Momentum Exchange |Δp|",
              title="RK4 + Courant: Momentum Exchange vs Impact Parameter",
              legend=:topright, grid=true, color=:blue)
    
    plot!(p1, b_values, dp_numerical,
          marker=:square, markersize=6, linewidth=2,
          label="Numerical (RK4)",
          linestyle=:dash, color=:red)
    
    # Add 1/b reference line
    b_ref = [minimum(b_values), maximum(b_values)]
    dp_ref = dp_analytical[1] * b_values[1] ./ b_ref  # Scale to match first point
    plot!(p1, b_ref, dp_ref,
          linewidth=1, linestyle=:dot, color=:gray,
          label="1/b scaling")
    
    # Create relative error plot
    p2 = plot(b_values, rel_errors,
              xscale=:log10, yscale=:log10,
              marker=:diamond, markersize=6, linewidth=2,
              color=:red,
              xlabel="Impact Parameter b", ylabel="Relative Error (%)",
              title="RK4 Numerical Accuracy vs Impact Parameter",
              legend=false, grid=true)
    
    # Combine plots
    combined_plot = plot(p1, p2, layout=(2,1), size=(800, 600))
    
    # Save plot
    savefig(combined_plot, "momentum_exchange_sweep_RK4.png")
    println("\n=== Plot Generated ===")
    println("Saved plot as 'momentum_exchange_sweep_RK4.png'")
    println("- Top panel: Momentum exchange scaling (shows perfect 1/b behavior)")
    println("- Bottom panel: Numerical accuracy across impact parameter range")
    println("- Uses RK4 integrator with adaptive timestepping (CFL condition)")
    
    return combined_plot
end

    # Run the momentum exchange sweep
    results = main()
    
    # Generate plots
    plot_momentum_exchange_rk4(results)
