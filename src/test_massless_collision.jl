using LinearAlgebra
using CSV
using ForwardDiff
using OrdinaryDiffEq
using Logging
using DoubleFloats
using Statistics

include("pomin.jl")
using .pomin

"""
    generate_massless_collision_data()

Generate initial data for two massless particles approaching each other head-on.
The particles are set up in their center-of-momentum frame with equal and opposite momenta.

Returns:
- Particles: struct containing masses, positions, and momenta
"""
function generate_massless_collision_data()
    # Masses (both massless)
    m1 = Double64(0.0)
    m2 = Double64(0.0)
    
    # Initial positions
    # Particle 1 on left, Particle 2 on right
    separation = Double64(100.0)  # Initial separation
    q1 = Double64[-separation/2, 0.0, 0.0]
    q2 = Double64[separation/2, 0.0, 0.0]
    
    # Initial momenta
    # Equal and opposite momenta for center-of-momentum frame
    p_mag = Double64(1.0)  # Magnitude of momentum
    p1 = Double64[p_mag, 0.0, 0.0]  # Moving right
    p2 = Double64[-p_mag, 0.0, 0.0] # Moving left
    
    return pomin.Particles([m1, m2], [q1, q2], [p1, p2])
end

"""
    massless_collision_eom!(du, u, params, t)

Equations of motion for two massless particles using the post-Minkowski Hamiltonian.
This is the form required by Julia's ODE solvers.
"""
function massless_collision_eom!(du, u, params, t)
    # Extract masses from params
    m1, m2 = params
    
    # Calculate gradient of Hamiltonian using ForwardDiff
    dH_val = pomin.dH(3, [m1, m2], u)
    
    # Convert to symplectic form
    du .= pomin.Jsympl(dH_val)
end
    
"""
    integrate_massless_collision(initial_data::pomin.Particles, tspan::Tuple{Real,Real}; 
                               abstol::Real=1e-12, reltol::Real=1e-12)

Integrate the massless particle collision using the PoMiN integrator interface.

Parameters:
- initial_data: Particles struct containing initial masses, positions, and momenta
- tspan: tuple of (start_time, end_time)
- abstol: absolute tolerance for the integrator
- reltol: relative tolerance for the integrator

Returns:
- ODESolution: Solution object containing the integrated trajectory
"""
function integrate_massless_collision(initial_data::pomin.Particles, tspan::Tuple{Real,Real}; 
                                    abstol::Real=1e-12, reltol::Real=1e-12)
    # Convert Particles to phase space vector
    z0 = vcat(vcat(initial_data.q...), vcat(initial_data.p...))
    
    # Convert tolerances to match z0 element type
    abstol_t = convert(eltype(z0), abstol)
    reltol_t = convert(eltype(z0), reltol)
    
    # Use the integrator interface
    return pomin.jlintegratorfull(z0, massless_collision_eom!, tspan, initial_data.m, abstol_t, reltol_t)
end

"""
    analyze_collision(solution::ODESolution)

Analyze the collision properties.
Returns a dictionary with various collision parameters.
"""
function analyze_collision(solution::ODESolution)
    # Get time and state vectors
    t = solution.t
    u = solution.u
    n_steps = length(t)
    
    # Extract positions and momenta for each timestep
    q1 = [u[i][1:3] for i in 1:n_steps]
    q2 = [u[i][4:6] for i in 1:n_steps]
    p1 = [u[i][7:9] for i in 1:n_steps]
    p2 = [u[i][10:12] for i in 1:n_steps]
    
    # Calculate separation at each timestep
    r = [q2[i] - q1[i] for i in 1:n_steps]
    r_mag = [norm(r[i]) for i in 1:n_steps]
    
    # Calculate total momentum
    p_tot = [p1[i] + p2[i] for i in 1:n_steps]
    p_tot_mag = [norm(p_tot[i]) for i in 1:n_steps]
    
    # Calculate relative momentum
    p_rel = [p2[i] - p1[i] for i in 1:n_steps]
    p_rel_mag = [norm(p_rel[i]) for i in 1:n_steps]
    
    # Find minimum separation
    min_sep = minimum(r_mag)
    min_sep_idx = argmin(r_mag)
    t_min_sep = t[min_sep_idx]
    
    return Dict(
        "Minimum separation" => min_sep,
        "Time of minimum separation" => t_min_sep,
        "Total momentum variation" => maximum(p_tot_mag),
        "Relative momentum variation" => (maximum(p_rel_mag) - minimum(p_rel_mag)) / mean(p_rel_mag),
        "Initial separation" => r_mag[1],
        "Final separation" => r_mag[end]
    )
end

"""
    print_analysis(collision_analysis)

Helper function to print collision analysis results in a consistent format.
"""
function print_analysis(collision_analysis)
    println("Collision Analysis:")
    println("  Initial separation: ", collision_analysis["Initial separation"])
    println("  Minimum separation: ", collision_analysis["Minimum separation"])
    println("  Time of minimum separation: ", collision_analysis["Time of minimum separation"])
    println("  Final separation: ", collision_analysis["Final separation"])
    println("\nConservation checks:")
    println("  Total momentum variation: ", collision_analysis["Total momentum variation"])
    println("  Relative momentum variation: ", collision_analysis["Relative momentum variation"])
end

# Main execution code
function main()
    # Generate initial data
    println("Generating initial conditions...")
    initial_data = generate_massless_collision_data()
    
    # Set integration parameters
    tstart = Double64(0.0)
    tend = Double64(200.0)  # Long enough to see the full collision
    
    # Integrate the system
    println("\nStarting integration...")
    sol = integrate_massless_collision(initial_data, (tstart, tend))
    
    println("\nIntegration complete:")
    println("  Number of timesteps: ", length(sol.t))
    println("  Time range: ", sol.t[1], " to ", sol.t[end])
    
    # Analyze collision
    println("\nAnalyzing collision properties...")
    analysis = analyze_collision(sol)
    print_analysis(analysis)
    
    return sol, analysis
end

# Run the simulation
println("Starting massless particle collision simulation...")
sol, analysis = main()
