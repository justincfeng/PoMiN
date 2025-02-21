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
    generate_sun_earth_initial_data()

Generate initial data for Sun-Earth system in G=c=1 units.
Uses parameters from the PoMiN test suite, which are already in geometrized units where G=c=1.

Returns:
- Particles: struct containing masses, positions, and momenta
"""
function generate_sun_earth_initial_data()
    # Masses in geometrized units (G=c=1)
    m_sun = Double64(1.0)  # Sun mass is unit mass
    m_earth = Double64(3.0033693739E-6)  # Earth mass in solar masses
    
    # Initial positions in geometrized units
    # Sun at origin
    q_sun = Double64[0.0, 0.0, 0.0]
    # Earth position (from test files)
    q_earth = Double64[99615330.0, 0.0, 0.0]
    
    # Initial momenta (from test files)
    p_sun = Double64[0.0, 0.0, 0.0]
    p_earth = Double64[0.0, 3.03440999999999E-10, 0.0]
    
    return pomin.Particles([m_sun, m_earth], [q_sun, q_earth], [p_sun, p_earth])
end

"""
    binary_system_eom!(du, u, params, t)

Equations of motion for the binary system using the post-Minkowski Hamiltonian.
This is the form required by Julia's ODE solvers.
"""
function binary_system_eom!(du, u, params, t)
    # Extract masses from params
    m_sun, m_earth = params
    
    # Calculate gradient of Hamiltonian
    dH_val = pomin.dH(3, [m_sun, m_earth], u)
    
    # Convert to symplectic form
    du .= pomin.Jsympl(dH_val)
end

"""
    integrate_binary_system(initial_data::pomin.Particles, G::Real, tspan::Tuple{Real,Real}; 
                          abstol::Real=1e-10, reltol::Real=1e-10)

Integrate a binary system using the PoMiN integrator interface.

Parameters:
- initial_data: Particles struct containing initial masses, positions, and momenta
- G: gravitational constant
- tspan: tuple of (start_time, end_time)
- abstol: absolute tolerance for the integrator
- reltol: relative tolerance for the integrator

Returns:
- ODESolution: Solution object containing the integrated trajectory
"""
function integrate_binary_system(initial_data::pomin.Particles, G::Real, tspan::Tuple{Real,Real}; 
                               abstol::Real=1e-10, reltol::Real=1e-10)
    # Convert Particles to phase space vector
    z0 = vcat(vcat(initial_data.q...), vcat(initial_data.p...))
    
    # Convert tolerances to match z0 element type
    abstol_t = convert(eltype(z0), abstol)
    reltol_t = convert(eltype(z0), reltol)
    
    # Use the new integrator interface
    return pomin.jlintegratorfull(z0, binary_system_eom!, tspan, initial_data.m, abstol_t, reltol_t)
end

"""
    integrate_binary_system_rk4(initial_data::pomin.Particles, G::Real, tspan::Tuple{Real,Real}; 
                              orbits_per_step::Real=10000)

Integrate a binary system using the RK4 integrator.

Parameters:
- initial_data: Particles struct containing initial masses, positions, and momenta
- G: gravitational constant
- tspan: tuple of (start_time, end_time)
- orbits_per_step: number of timesteps per expected orbit (default: 10000)

Returns:
- soln: solution struct containing the integrated trajectory
"""
function integrate_binary_system_rk4(initial_data::pomin.Particles, G::Real, tspan::Tuple{Real,Real}; 
                                   orbits_per_step::Real=10000)
    # Convert Particles to phase space vector
    z0 = vcat(vcat(initial_data.q...), vcat(initial_data.p...))
    
    # Create function that returns dH/dz
    function dz(z)
        return pomin.dH(3, initial_data.m, z)
    end
    
    # Time adaptation function - just returns constant timestep
    function tadapt(dt, z, dz)
        return convert(eltype(z), dt)
    end
    
    # Estimate orbital period using Kepler's third law
    # a ≈ initial radius for nearly circular orbit
    r0 = norm(initial_data.q[2] - initial_data.q[1])
    M = initial_data.m[1] + initial_data.m[2]  # Total mass
    T = 2π * sqrt(r0^3 / (G * M))  # Full Kepler period formula
    
    # Set timestep based on estimated period
    dt = T / orbits_per_step
    
    # Number of timesteps
    n_steps = Int(ceil((tspan[2] - tspan[1])/dt))
    
    println("  Estimated orbital period: ", T)
    println("  Timestep: ", dt)
    println("  Steps per orbit: ", orbits_per_step)
    println("  Total steps: ", n_steps)
    
    # Use RK4 integrator
    return pomin.hrkintegrator(3, 2, z0, dz, dt, tadapt, tspan, n_steps)
end

"""
    integrate_binary_system_julia(initial_data::pomin.Particles, G::Real, tspan::Tuple{Real,Real}; 
                                abstol::Real=1e-10, reltol::Real=1e-10)

Integrate a binary system using Julia's ODE solvers via the direct interface.

Parameters:
- initial_data: Particles struct containing initial masses, positions, and momenta
- G: gravitational constant
- tspan: tuple of (start_time, end_time)
- abstol: absolute tolerance for the integrator
- reltol: relative tolerance for the integrator

Returns:
- ODESolution: Solution object containing the integrated trajectory
"""
function integrate_binary_system_julia(initial_data::pomin.Particles, G::Real, tspan::Tuple{Real,Real}; 
                                     abstol::Real=1e-10, reltol::Real=1e-10)
    # Convert Particles to phase space vector
    z0 = vcat(vcat(initial_data.q...), vcat(initial_data.p...))
    
    # Convert tolerances to match z0 element type
    abstol_t = convert(eltype(z0), abstol)
    reltol_t = convert(eltype(z0), reltol)
    
    # Create and solve ODE problem directly
    prob = ODEProblem(binary_system_eom!, z0, tspan, initial_data.m)
    return solve(prob, Tsit5(), abstol=abstol_t, reltol=reltol_t)
end

"""
    analyze_orbit(solution::Union{ODESolution,pomin.soln}, m_sun::T, m_earth::T) where T<:Real

Analyze the orbital properties to verify Keplerian nature.
Returns a dictionary with various orbital parameters.
"""
function analyze_orbit(solution::Union{ODESolution,pomin.soln}, m_sun::T, m_earth::T) where T<:Real
    # Get time and state vectors based on solution type
    if solution isa ODESolution
        t = solution.t
        u = solution.u
    else
        t = solution.t
        u = solution.z
    end
    
    n_steps = length(t)
    
    # Extract positions and momenta from solution
    # Each u[i] is the full state vector at time t[i]
    q_sun = [[u[i][j] for j in 1:3] for i in 1:n_steps]
    q_earth = [[u[i][j] for j in 4:6] for i in 1:n_steps]
    p_sun = [[u[i][j] for j in 7:9] for i in 1:n_steps]
    p_earth = [[u[i][j] for j in 10:12] for i in 1:n_steps]
    
    # Calculate relative position and velocity at each timestep
    r = [q_earth[i] - q_sun[i] for i in 1:n_steps]
    v = [p_earth[i]/m_earth - p_sun[i]/m_sun for i in 1:n_steps]
    
    # Calculate angular momentum per unit mass (should be conserved)
    L = [cross(r[i], v[i]) for i in 1:n_steps]
    L_mag = [norm(L[i]) for i in 1:n_steps]
    L_var = (maximum(L_mag) - minimum(L_mag)) / mean(L_mag)
    
    # Calculate energy per unit mass (should be conserved)
    # E = v²/2 - GM/r
    E = [dot(v[i],v[i])/2 - 1/norm(r[i]) for i in 1:n_steps]
    E_var = (maximum(E) - minimum(E)) / abs(mean(E))
    
    # Calculate eccentricity vector (should be constant for Keplerian orbit)
    # e = v × L - r̂
    e = [cross(v[i], L[i]) / 1.0 - r[i]/norm(r[i]) for i in 1:n_steps]
    e_mag = [norm(e[i]) for i in 1:n_steps]
    e_var = (maximum(e_mag) - minimum(e_mag)) / mean(e_mag)
    
    # Find orbital period by looking for returns to initial position
    r_mag = [norm(r[i]) for i in 1:n_steps]
    
    # Find peaks in r_mag to estimate period
    periods = T[]
    last_peak_idx = 1
    for i in 2:n_steps-1
        if r_mag[i] > r_mag[i-1] && r_mag[i] > r_mag[i+1]
            if last_peak_idx > 1
                # Use actual time difference instead of assuming constant dt
                push!(periods, t[i] - t[last_peak_idx])
            end
            last_peak_idx = i
        end
    end
    
    # Calculate semi-major axis from energy
    # Use Complex sqrt to handle negative energies
    a = abs(-1/(2*mean(E)))
    
    # Calculate expected period from Kepler's third law
    # T² = 4π²a³/GM with G=1, M≈m_sun (since m_earth << m_sun)
    expected_period = 2π * sqrt(Complex(a^3))
    
    return Dict(
        "Angular momentum variation" => L_var,
        "Energy variation" => E_var,
        "Eccentricity variation" => e_var,
        "Mean eccentricity" => mean(e_mag),
        "Semi-major axis" => a,
        "Orbital periods" => periods,
        "Expected period" => abs(expected_period),  # Take absolute value of complex result
        "Number of orbits" => length(periods)
    )
end

"""
    print_analysis(orbital_analysis)

Helper function to print orbital analysis results in a consistent format.
"""
function print_analysis(orbital_analysis)
    println("Relative variations (should be small for Keplerian orbit):")
    println("  Angular momentum: ", orbital_analysis["Angular momentum variation"])
    println("  Energy: ", orbital_analysis["Energy variation"])
    println("  Eccentricity: ", orbital_analysis["Eccentricity variation"])
    println("\nOrbital Parameters:")
    println("  Mean eccentricity: ", orbital_analysis["Mean eccentricity"])
    println("  Semi-major axis: ", orbital_analysis["Semi-major axis"])
    println("  Number of completed orbits: ", orbital_analysis["Number of orbits"])
    if !isempty(orbital_analysis["Orbital periods"])
        println("  Mean orbital period: ", mean(orbital_analysis["Orbital periods"]))
        println("  Expected period (Kepler's Third Law): ", orbital_analysis["Expected period"])
        println("  Period variation: ", (maximum(orbital_analysis["Orbital periods"]) - minimum(orbital_analysis["Orbital periods"])) / mean(orbital_analysis["Orbital periods"]))
    else
        println("  Not enough orbits completed to calculate period")
    end
end

# Main execution code
function main()
    # Generate initial data using Sun-Earth parameters
    println("Generating initial conditions...")
    initial_data = generate_sun_earth_initial_data()
    
    # Set integration parameters
    tstart = Double64(0.0)
    tend = Double64(6.4E+13)  # About 10 Earth years
    
    # Integrate the system using different methods
    println("\nStarting integrations...")
    
    println("\n1. Using PoMiN integrator interface...")
    sol_pomin = integrate_binary_system(initial_data, Double64(1.0), (tstart, tend))
    
    println("\n2. Using RK4 integrator...")
    sol_rk4 = integrate_binary_system_rk4(initial_data, Double64(1.0), (tstart, tend))
    
    println("\n3. Using direct Julia integrator...")
    sol_julia = integrate_binary_system_julia(initial_data, Double64(1.0), (tstart, tend))
    
    # Print comparison of methods
    println("\nComparison of integration methods:")
    println("1. PoMiN interface:")
    println("   Number of timesteps: ", length(sol_pomin.t))
    println("   Time range: ", sol_pomin.t[1], " to ", sol_pomin.t[end])
    
    println("\n2. RK4:")
    println("   Number of timesteps: ", length(sol_rk4.t))
    println("   Time range: ", sol_rk4.t[1], " to ", sol_rk4.t[end])
    
    println("\n3. Direct Julia:")
    println("   Number of timesteps: ", length(sol_julia.t))
    println("   Time range: ", sol_julia.t[1], " to ", sol_julia.t[end])
    
    # Analyze orbits
    println("\nAnalyzing orbital properties...")
    
    println("\n1. PoMiN interface:")
    analysis_pomin = analyze_orbit(sol_pomin, Double64(1.0), Double64(3.0033693739E-6))
    print_analysis(analysis_pomin)
    
    println("\n2. RK4:")
    analysis_rk4 = analyze_orbit(sol_rk4, Double64(1.0), Double64(3.0033693739E-6))
    print_analysis(analysis_rk4)
    
    println("\n3. Direct Julia:")
    analysis_julia = analyze_orbit(sol_julia, Double64(1.0), Double64(3.0033693739E-6))
    print_analysis(analysis_julia)
    
    return (sol_pomin, sol_rk4, sol_julia), 
           (analysis_pomin, analysis_rk4, analysis_julia)
end

# Run the simulation
println("Starting Sun-Earth simulation comparison...")
sols, analyses = main()
