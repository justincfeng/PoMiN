include("../pomin.jl")
using LinearAlgebra
using Plots

"""
    dpfunc_analytical(p, b, m1, m2; G=1, c=1)

Analytical momentum exchange formula for post-Minkowskian scattering.
Returns the change in momentum for particle 1 in the y-direction.
"""
function dpfunc_analytical(p, b, m1, m2; G=1, c=1)
    E1 = c*√((c*m1)^2 + p^2)
    E2 = c*√((c*m2)^2 + p^2)
    dp = (2*G*(E1*E2)^2)/(b*p*(E1+E2)) *
         (1 + (1/E1^2 + 1/E2^2 + 4/(E1*E2))*p^2 + p^4/(E1*E2)^2)
    return dp
end

"""
    momentum_exchange_sweep()

Perform momentum exchange calculations over a wide range of impact parameters.
Uses the efficient Nrec=-1 option to save only initial and final states.
"""
function momentum_exchange_sweep()
    println("=== PoMiN Momentum Exchange Sweep ===")
    println("Testing impact parameter scaling with Nrec=-1 (save-last-only)\n")
    
    # Physical parameters
    m1 = 1.0
    m2 = 1.0  
    p = 10.0
    
    # Impact parameter sweep parameters
    b_base = 10.0
    scale_factor = 1.5  # Geometric progression factor
    nn_values = [5, 10, 15, 20, 25]  # Range of scaling exponents
    
    results = []
    
    for nn in nn_values
        println("Testing with nn = $nn")
        
        # Generate impact parameter with wide dynamic range
        b = b_base * scale_factor^nn
        dx = 100*b  # Initial separation
        
        println("  Impact parameter b = $(round(b, sigdigits=6))")
        println("  Initial separation dx = $(round(dx, sigdigits=6))")
        
        # Set up scattering system
        system = pomin.setup_scattering(m1, m2, p, b, dx)
        
        # Solve with Julia ODE integrator (save only last point)
        sol = pomin.solve(system, pomin.ParametersJulia((0.0, 2*dx), 
                         integrator="Tsit5", atol=1e-12, rtol=1e-12, Nrec=-1))
        
        # Extract momentum change (particle 1, y-component)
        # Initial momentum: sol.u[1] = [q1x, q1y, q1z, q2x, q2y, q2z, p1x, p1y, p1z, p2x, p2y, p2z]
        # Final momentum:   sol.u[end]
        p1y_initial = sol.u[1][8]   # p1y component (8th element)
        p1y_final = sol.u[end][8]   # p1y component at end
        dp_numerical = p1y_final - p1y_initial
        
        # Analytical prediction
        dp_analytical = dpfunc_analytical(p, b, m1, m2)
        
        # Calculate relative error
        rel_error = abs(dp_numerical - dp_analytical) / abs(dp_analytical) * 100
        
        println("  Analytical Δp = $(round(dp_analytical, sigdigits=8))")
        println("  Numerical  Δp = $(round(dp_numerical, sigdigits=8))")
        println("  Relative error = $(round(rel_error, digits=3))%")
        println("  Integration time = $(round(sol.t[end], sigdigits=6))")
        println("  Saved $(length(sol.t)) time points (memory efficient!)\n")
        
        # Store results
        push!(results, (nn=nn, b=b, dp_analytical=dp_analytical, 
                       dp_numerical=dp_numerical, rel_error=rel_error,
                       final_time=sol.t[end]))
    end
    
    # Summary
    println("=== Summary ===")
    println("nn\tb\t\tΔp_analytical\tΔp_numerical\tRel_Error(%)")
    println("-"^70)
    for r in results
        println("$(r.nn)\t$(round(r.b, sigdigits=4))\t\t$(round(r.dp_analytical, sigdigits=6))\t$(round(r.dp_numerical, sigdigits=6))\t$(round(r.rel_error, digits=2))")
    end
    
    return results
end

"""
    plot_momentum_exchange(results)

Generate plots showing momentum exchange scaling and accuracy.
"""
function plot_momentum_exchange(results)
    # Extract data for plotting
    b_values = [r.b for r in results]
    dp_analytical = [r.dp_analytical for r in results]
    dp_numerical = [r.dp_numerical for r in results]
    rel_errors = [r.rel_error for r in results]
    
    # Create log-log plot of momentum exchange vs impact parameter
    p1 = plot(b_values, dp_analytical, 
              xscale=:log10, yscale=:log10,
              marker=:circle, markersize=6, linewidth=2,
              label="Analytical (PM Theory)",
              xlabel="Impact Parameter b", ylabel="Momentum Exchange |Δp|",
              title="Momentum Exchange vs Impact Parameter",
              legend=:topright, grid=true)
    
    plot!(p1, b_values, dp_numerical,
          marker=:square, markersize=6, linewidth=2,
          label="Numerical (PoMiN)",
          linestyle=:dash)
    
    # Add 1/b reference line
    b_ref = [minimum(b_values), maximum(b_values)]
    dp_ref = dp_analytical[1] * b_values[1] ./ b_ref  # Scale to match first point
    plot!(p1, b_ref, dp_ref,
          linewidth=1, linestyle=:dot, color=:gray,
          label="1/b scaling")
    
    # Create relative error plot
    p2 = plot(b_values, rel_errors,
              xscale=:log10,
              marker=:diamond, markersize=6, linewidth=2,
              color=:red,
              xlabel="Impact Parameter b", ylabel="Relative Error (%)",
              title="Numerical Accuracy vs Impact Parameter",
              legend=false, grid=true)
    
    # Combine plots
    combined_plot = plot(p1, p2, layout=(2,1), size=(800, 600))
    
    # Save plot
    savefig(combined_plot, "momentum_exchange_sweep.png")
    println("\n=== Plot Generated ===")
    println("Saved plot as 'momentum_exchange_sweep.png'")
    println("- Top panel: Momentum exchange scaling (shows perfect 1/b behavior)")
    println("- Bottom panel: Numerical accuracy across impact parameter range")
    
    return combined_plot
end

# Run the momentum exchange sweep
results = momentum_exchange_sweep()

# Generate plots
plot_momentum_exchange(results)

