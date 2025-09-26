include("../pomin.jl")
using LinearAlgebra
using Plots
using LaTeXStrings
using DoubleFloats

#-----------------------------------------------------------------------
"""
    momentum_exchange_sweep()

Perform momentum exchange calculations over a wide range of impact
parameters. Uses the efficient Nrec=-1 option to save only initial and
final states.
"""
function momentum_exchange_sweep()
   
    # Physical parameters (using DoubleFloats for higher precision)
    m1 = Double64(1.0)  # Massless particle 1
    m2 = Double64(0.0)  # Massless particle 2
    p = Double64(10.0)
    
    # Impact parameter sweep parameters
    b_base = Double64(10.0)
    scale_factor = Double64("1.1108305558745590335689712446765042841434478759765625")
    nn_values = 50:20:350  # Range from 1 to 250 in steps of 50
    
    results = []
    
    for nn in nn_values       
        # Generate impact parameter with wide dynamic range
        b = b_base * scale_factor^nn
        dx = Double64(1e9)*b

        # Set up scattering system (using DoubleFloats type)
        system = pomin.to_com_frame(pomin.setup_scattering(m1, m2, p, b, 
                                                     dx, tpfl=Double64))
        
        # Calculate duration using alternative MXIC strategy
        # For massless particles, velocity = c = 1 (natural units)
        v1 = p / sqrt(m1^2 + p^2)  # Speed for massive particle
        v2 = Double64(1.0)         # Speed of light, massless particle
        τ = dx / (v1 + v2)         # Time for particles to meet
        t_flight = Double64(10.0) * τ         # Total scattering time
        # Solve with Julia ODE integrator (save only last point)
        sol = pomin.solve(system, 
                          pomin.ParametersJulia(
                            (Double64(0.0),t_flight), 
                            integrator="Vern9", 
                            atol=1e-18, rtol=1e-18, Nrec=-1))
        # Extract momentum change (particle 1, y-component)
        # Initial momentum: sol.u[1] = [q1x,q1y,q1z,q2x,q2y,q2z,
        #                               p1x,p1y,p1z,p2x,p2y,p2z]
        # Final momentum:   sol.u[end]
        p1y_initial = sol.u[1][8]   # p1y component (8th element)
        p1y_final = sol.u[end][8]   # p1y component at end
        dp_numerical = p1y_final - p1y_initial
        
        # Analytical prediction
        dp_analytical = pomin.HamPM.dp_scatter(p, b, m1, m2)
        
        # Calculate relative error
        rel_error = abs(abs(dp_numerical) - abs(dp_analytical)) / 
                        abs(dp_analytical) * 100
        
        # Store results
        push!(results, (nn=nn, b=b, dp_analytical=dp_analytical, 
                       dp_numerical=dp_numerical, rel_error=rel_error,
                       final_time=sol.t[end]))
    end

    return results
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    plot_momentum_exchange_scaling(results)

Generate publication-quality momentum exchange scaling plot for Physical Review D.
"""
function plot_momentum_exchange_scaling(results)
    # Extract data for plotting
    b_values = [r.b for r in results]
    dp_analytical = [r.dp_analytical for r in results]
    dp_numerical = [r.dp_numerical for r in results]
    
    # Set publication-quality defaults
    default(fontfamily="Computer Modern", 
            guidefontsize=10, tickfontsize=8, legendfontsize=9,
            linewidth=2, markersize=4, markerstrokewidth=0,
            dpi=300, size=(400, 300))
    
    # Create log-log plot of momentum exchange vs impact parameter
    p1 = plot(b_values, abs.(dp_analytical), 
              xscale=:log10, yscale=:log10,
              marker=:circle, markersize=4, linewidth=2,
              color=:black, markerstrokewidth=0,
              label="Analytical",
              xlabel=L"Impact parameter $b$", 
              ylabel=L"Momentum exchange $|\Delta p|$",
              legend=:topright, grid=true, minorgrid=false,
              framestyle=:box,
              xticks=10.0 .^ (1:2:60), yticks=10.0 .^ (-20:2:-4))
    
    plot!(p1, b_values, abs.(dp_numerical),
          marker=:circle, markersize=4, linewidth=2,
          color=:red, markerstrokewidth=0,
          label="Numerical",
          linestyle=:solid)
    
    # Add 1/b reference line
    b_ref = [minimum(b_values), maximum(b_values)]
    dp_ref = abs(dp_analytical[1]) * b_values[1] ./ b_ref
    plot!(p1, b_ref, dp_ref,
          linewidth=1.5, linestyle=:dash, color=:gray,
          label=L"$\propto b^{-1}$")

    return plot(p1, left_margin=5Plots.mm, bottom_margin=5Plots.mm,
                top_margin=2Plots.mm, right_margin=2Plots.mm)
end

"""
    plot_momentum_exchange_error(results)

Generate publication-quality relative error plot for Physical Review D.
"""
function plot_momentum_exchange_error(results)
    # Extract data for plotting
    b_values = [r.b for r in results]
    rel_errors = [r.rel_error for r in results]
    
    # Set publication-quality defaults
    default(fontfamily="Computer Modern", 
            guidefontsize=10, tickfontsize=8, legendfontsize=9,
            linewidth=2, markersize=4, markerstrokewidth=0,
            dpi=300, size=(400, 300))
    
    # Create relative error plot
    p2 = plot(b_values, rel_errors,
              xscale=:log10, yscale=:log10,
              marker=:circle, markersize=4, linewidth=2,
              color=:steelblue, markerstrokewidth=0,
              xlabel=L"Impact parameter $b$", 
              ylabel="Relative error (%)",
              legend=false, grid=true, minorgrid=false,
              framestyle=:box,
              xticks=10.0 .^ (1:2:60), yticks=10.0 .^ (-12:2:2))

    return plot(p2, left_margin=5Plots.mm, bottom_margin=5Plots.mm,
                top_margin=2Plots.mm, right_margin=2Plots.mm)
end

# Run the momentum exchange sweep
results = momentum_exchange_sweep()

# Generate separate plots
scaling_plot = plot_momentum_exchange_scaling(results)
error_plot = plot_momentum_exchange_error(results)

# Save separate plots
savefig(scaling_plot, "momentum_exchange_scaling_mixed.pdf")
savefig(error_plot, "momentum_exchange_error_mixed.pdf")
